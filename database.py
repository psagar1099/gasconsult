"""
Database module for persistent chat history storage.

This module provides SQLite-based storage for conversations and messages.
It's designed to be completely isolated from the main app - if database
operations fail, the app continues working with session-based storage.

Production database location: /var/lib/gasconsult/gasconsult.db
Development fallback: ./gasconsult.db (current directory)
"""

import sqlite3
import json
import uuid
import os
import logging
from datetime import datetime
from contextlib import contextmanager
from typing import List, Dict, Optional, Any

# Configure logging
logger = logging.getLogger(__name__)

# Database configuration
DB_DIR = os.getenv('GASCONSULT_DB_DIR', '/var/lib/gasconsult')
DB_NAME = 'gasconsult.db'
DB_PATH = os.path.join(DB_DIR, DB_NAME)

# Fallback to current directory if production path doesn't exist
# (for development/testing)
if not os.path.exists(DB_DIR):
    logger.warning(f"⚠️  Production DB directory {DB_DIR} doesn't exist!")
    logger.warning(f"⚠️  Using fallback path. Data will NOT persist across redeployments!")
    logger.warning(f"⚠️  To fix: Add a persistent disk in Render mounted at {DB_DIR}")
    DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)
    logger.info(f"Using fallback database path: {DB_PATH}")
else:
    logger.info(f"✓ Using persistent database at: {DB_PATH}")


@contextmanager
def get_db_connection():
    """
    Context manager for database connections.
    Ensures connections are properly closed even if errors occur.

    Usage:
        with get_db_connection() as conn:
            conn.execute(...)
    """
    conn = None
    try:
        conn = sqlite3.connect(DB_PATH, timeout=10.0)
        conn.row_factory = sqlite3.Row  # Enable column access by name
        yield conn
        conn.commit()
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Database error: {e}")
        raise
    finally:
        if conn:
            conn.close()


def init_db():
    """
    Initialize the database with required tables and indexes.
    Safe to call multiple times - uses IF NOT EXISTS.

    Returns:
        bool: True if successful, False if error occurred
    """
    try:
        # Ensure directory exists (for production path)
        db_dir = os.path.dirname(DB_PATH)
        if db_dir and not os.path.exists(db_dir):
            try:
                os.makedirs(db_dir, mode=0o755)
                logger.info(f"Created database directory: {db_dir}")
            except PermissionError:
                logger.error(f"Permission denied creating {db_dir}, using fallback")
                # Will use fallback path set at module level
                return False

        with get_db_connection() as conn:
            # Create users table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS users (
                    id TEXT PRIMARY KEY,
                    email TEXT UNIQUE NOT NULL,
                    password_hash TEXT,
                    full_name TEXT,
                    oauth_provider TEXT CHECK(oauth_provider IN ('google', 'apple', NULL)),
                    oauth_provider_id TEXT,
                    is_verified BOOLEAN DEFAULT 0,
                    verification_token TEXT,
                    reset_token TEXT,
                    reset_token_expires TIMESTAMP,
                    subscription_tier TEXT DEFAULT 'free' CHECK(subscription_tier IN ('free', 'pro', 'team')),
                    subscription_status TEXT DEFAULT 'active' CHECK(subscription_status IN ('active', 'cancelled', 'expired')),
                    stripe_customer_id TEXT,
                    stripe_subscription_id TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    last_login TIMESTAMP,
                    is_active BOOLEAN DEFAULT 1,
                    UNIQUE(oauth_provider, oauth_provider_id)
                )
            """)

            # Create conversations table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS conversations (
                    id TEXT PRIMARY KEY,
                    user_session_id TEXT,
                    user_id TEXT,
                    title TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    message_count INTEGER DEFAULT 0,
                    is_active BOOLEAN DEFAULT 1,
                    FOREIGN KEY (user_id) REFERENCES users(id) ON DELETE CASCADE
                )
            """)

            # Create messages table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS messages (
                    id TEXT PRIMARY KEY,
                    conversation_id TEXT NOT NULL,
                    role TEXT NOT NULL CHECK(role IN ('user', 'assistant')),
                    content TEXT NOT NULL,
                    paper_references TEXT,
                    num_papers INTEGER DEFAULT 0,
                    evidence_strength TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (conversation_id) REFERENCES conversations(id) ON DELETE CASCADE
                )
            """)

            # Create indexes for performance
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_users_email
                ON users(email)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_users_verification_token
                ON users(verification_token)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_users_reset_token
                ON users(reset_token)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_conversations_session
                ON conversations(user_session_id, created_at DESC)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_conversations_user
                ON conversations(user_id, created_at DESC)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_conversations_active
                ON conversations(user_session_id, is_active, created_at DESC)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_messages_conversation
                ON messages(conversation_id, created_at ASC)
            """)

            # Create bookmarks table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS bookmarks (
                    id TEXT PRIMARY KEY,
                    user_session_id TEXT,
                    user_id TEXT,
                    query TEXT NOT NULL,
                    answer TEXT NOT NULL,
                    paper_references TEXT,
                    num_papers INTEGER DEFAULT 0,
                    evidence_strength TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (user_id) REFERENCES users(id) ON DELETE CASCADE
                )
            """)

            # Create shared_links table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS shared_links (
                    id TEXT PRIMARY KEY,
                    share_id TEXT UNIQUE NOT NULL,
                    query TEXT NOT NULL,
                    answer TEXT NOT NULL,
                    paper_references TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    expires_at TIMESTAMP NOT NULL
                )
            """)

            # Create indexes for bookmarks
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_bookmarks_user_session
                ON bookmarks(user_session_id, created_at DESC)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_bookmarks_user
                ON bookmarks(user_id, created_at DESC)
            """)

            # Create indexes for shared_links
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_shared_links_share_id
                ON shared_links(share_id)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_shared_links_expires
                ON shared_links(expires_at)
            """)

            logger.info(f"Database initialized successfully at {DB_PATH}")
            return True

    except Exception as e:
        logger.error(f"Failed to initialize database: {e}")
        return False


def save_conversation(
    id: str,
    user_session_id: str,
    title: str
) -> bool:
    """
    Create a new conversation record.

    Args:
        id: Unique conversation ID (UUID)
        user_session_id: Session ID or user ID
        title: Conversation title (auto-generated from first query)

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                INSERT INTO conversations (id, user_session_id, title, message_count)
                VALUES (?, ?, ?, 0)
            """, (id, user_session_id, title))

        logger.info(f"Saved conversation {id}: '{title}'")
        return True

    except Exception as e:
        logger.error(f"Failed to save conversation {id}: {e}")
        return False


def save_message(
    conversation_id: str,
    role: str,
    content: str,
    references: Optional[List[Dict]] = None,
    num_papers: int = 0,
    evidence_strength: Optional[str] = None
) -> bool:
    """
    Save a message to an existing conversation.

    Args:
        conversation_id: ID of the conversation
        role: 'user' or 'assistant'
        content: Message content
        references: List of paper references (for assistant messages)
        num_papers: Number of papers cited
        evidence_strength: 'high', 'moderate', or 'low'

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        message_id = str(uuid.uuid4())
        references_json = json.dumps(references) if references else None
        evidence_strength_json = json.dumps(evidence_strength) if evidence_strength else None

        with get_db_connection() as conn:
            # Insert message
            conn.execute("""
                INSERT INTO messages
                (id, conversation_id, role, content, paper_references, num_papers, evidence_strength)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (message_id, conversation_id, role, content, references_json, num_papers, evidence_strength_json))

            # Update conversation metadata
            conn.execute("""
                UPDATE conversations
                SET message_count = message_count + 1,
                    updated_at = CURRENT_TIMESTAMP
                WHERE id = ?
            """, (conversation_id,))

        logger.debug(f"Saved {role} message to conversation {conversation_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to save message to conversation {conversation_id}: {e}")
        return False


def get_conversations(
    user_session_id: str,
    limit: int = 50,
    offset: int = 0,
    include_inactive: bool = False
) -> List[Dict[str, Any]]:
    """
    Get list of conversations for a user.

    Args:
        user_session_id: Session ID or user ID
        limit: Maximum number of conversations to return
        offset: Offset for pagination
        include_inactive: Include soft-deleted conversations

    Returns:
        List of conversation dictionaries
    """
    try:
        with get_db_connection() as conn:
            query = """
                SELECT id, title, created_at, updated_at, message_count
                FROM conversations
                WHERE user_session_id = ?
            """

            if not include_inactive:
                query += " AND is_active = 1"

            query += " ORDER BY updated_at DESC LIMIT ? OFFSET ?"

            cursor = conn.execute(query, (user_session_id, limit, offset))
            rows = cursor.fetchall()

            conversations = []
            for row in rows:
                conversations.append({
                    'id': row['id'],
                    'title': row['title'],
                    'created_at': row['created_at'],
                    'updated_at': row['updated_at'],
                    'message_count': row['message_count']
                })

            return conversations

    except Exception as e:
        logger.error(f"Failed to get conversations for user {user_session_id}: {e}")
        return []


def get_conversation(conversation_id: str) -> Optional[Dict[str, Any]]:
    """
    Get a specific conversation with all its messages.

    Args:
        conversation_id: ID of the conversation

    Returns:
        Dictionary with conversation metadata and messages, or None if not found
    """
    try:
        with get_db_connection() as conn:
            # Get conversation metadata
            cursor = conn.execute("""
                SELECT id, user_session_id, title, created_at, updated_at, message_count
                FROM conversations
                WHERE id = ? AND is_active = 1
            """, (conversation_id,))

            conv_row = cursor.fetchone()
            if not conv_row:
                return None

            # Get all messages for this conversation
            cursor = conn.execute("""
                SELECT role, content, paper_references, num_papers, evidence_strength, created_at
                FROM messages
                WHERE conversation_id = ?
                ORDER BY created_at ASC
            """, (conversation_id,))

            messages = []
            for msg_row in cursor.fetchall():
                message = {
                    'role': msg_row['role'],
                    'content': msg_row['content']
                }

                # Add assistant-specific fields
                if msg_row['role'] == 'assistant':
                    message['references'] = json.loads(msg_row['paper_references']) if msg_row['paper_references'] else []
                    message['num_papers'] = msg_row['num_papers'] or 0
                    if msg_row['evidence_strength']:
                        message['evidence_strength'] = json.loads(msg_row['evidence_strength'])

                messages.append(message)

            return {
                'id': conv_row['id'],
                'user_session_id': conv_row['user_session_id'],
                'title': conv_row['title'],
                'created_at': conv_row['created_at'],
                'updated_at': conv_row['updated_at'],
                'message_count': conv_row['message_count'],
                'messages': messages
            }

    except Exception as e:
        logger.error(f"Failed to get conversation {conversation_id}: {e}")
        return None


def delete_conversation(conversation_id: str, hard_delete: bool = False) -> bool:
    """
    Delete a conversation (soft delete by default).

    Args:
        conversation_id: ID of the conversation to delete
        hard_delete: If True, permanently delete; if False, just mark inactive

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            if hard_delete:
                # Permanently delete (CASCADE will delete messages too)
                conn.execute("DELETE FROM conversations WHERE id = ?", (conversation_id,))
                logger.info(f"Hard deleted conversation {conversation_id}")
            else:
                # Soft delete (mark as inactive)
                conn.execute(
                    "UPDATE conversations SET is_active = 0 WHERE id = ?",
                    (conversation_id,)
                )
                logger.info(f"Soft deleted conversation {conversation_id}")

            return True

    except Exception as e:
        logger.error(f"Failed to delete conversation {conversation_id}: {e}")
        return False


def update_conversation_title(conversation_id: str, new_title: str) -> bool:
    """
    Update the title of a conversation.

    Args:
        conversation_id: ID of the conversation
        new_title: New title for the conversation

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE conversations
                SET title = ?, updated_at = CURRENT_TIMESTAMP
                WHERE id = ?
            """, (new_title, conversation_id))

        logger.info(f"Updated title for conversation {conversation_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to update conversation title {conversation_id}: {e}")
        return False


def search_conversations(
    user_session_id: str,
    search_query: str,
    limit: int = 20
) -> List[Dict[str, Any]]:
    """
    Search conversations by title or content.

    Args:
        user_session_id: Session ID or user ID
        search_query: Text to search for
        limit: Maximum number of results

    Returns:
        List of matching conversations
    """
    try:
        with get_db_connection() as conn:
            # Search in both conversation titles and message content
            cursor = conn.execute("""
                SELECT DISTINCT c.id, c.title, c.created_at, c.updated_at, c.message_count
                FROM conversations c
                LEFT JOIN messages m ON c.id = m.conversation_id
                WHERE c.user_session_id = ?
                  AND c.is_active = 1
                  AND (c.title LIKE ? OR m.content LIKE ?)
                ORDER BY c.updated_at DESC
                LIMIT ?
            """, (user_session_id, f'%{search_query}%', f'%{search_query}%', limit))

            conversations = []
            for row in cursor.fetchall():
                conversations.append({
                    'id': row['id'],
                    'title': row['title'],
                    'created_at': row['created_at'],
                    'updated_at': row['updated_at'],
                    'message_count': row['message_count']
                })

            return conversations

    except Exception as e:
        logger.error(f"Failed to search conversations: {e}")
        return []


def get_database_stats() -> Dict[str, Any]:
    """
    Get database statistics for monitoring.

    Returns:
        Dictionary with database stats
    """
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("SELECT COUNT(*) as total FROM conversations WHERE is_active = 1")
            total_conversations = cursor.fetchone()['total']

            cursor = conn.execute("SELECT COUNT(*) as total FROM messages")
            total_messages = cursor.fetchone()['total']

            cursor = conn.execute("""
                SELECT COUNT(DISTINCT user_session_id) as total
                FROM conversations
                WHERE is_active = 1
            """)
            total_users = cursor.fetchone()['total']

            # Get database file size
            db_size = os.path.getsize(DB_PATH) if os.path.exists(DB_PATH) else 0

            return {
                'database_path': DB_PATH,
                'database_size_mb': round(db_size / (1024 * 1024), 2),
                'total_conversations': total_conversations,
                'total_messages': total_messages,
                'total_users': total_users,
                'initialized': True
            }

    except Exception as e:
        logger.error(f"Failed to get database stats: {e}")
        return {
            'database_path': DB_PATH,
            'initialized': False,
            'error': str(e)
        }


# Utility function for generating conversation titles
def generate_conversation_title(first_query: str, max_length: int = 60) -> str:
    """
    Generate a conversation title from the first user query.

    Args:
        first_query: The first user message
        max_length: Maximum title length

    Returns:
        Generated title string
    """
    # Remove extra whitespace
    title = ' '.join(first_query.split())

    # Truncate if too long
    if len(title) > max_length:
        title = title[:max_length].rsplit(' ', 1)[0] + '...'

    return title


# ====== User Management Functions ======

def create_user(email: str, password_hash: str, full_name: Optional[str] = None, verification_token: Optional[str] = None) -> Optional[str]:
    """
    Create a new user account.

    Args:
        email: User's email address
        password_hash: Bcrypt hashed password
        full_name: User's full name (optional)
        verification_token: Email verification token (optional)

    Returns:
        User ID if successful, None otherwise
    """
    try:
        user_id = str(uuid.uuid4())
        with get_db_connection() as conn:
            conn.execute("""
                INSERT INTO users (id, email, password_hash, full_name, verification_token, is_verified)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (user_id, email.lower(), password_hash, full_name, verification_token, 1 if not verification_token else 0))

        logger.info(f"Created user account: {email}")
        return user_id

    except sqlite3.IntegrityError:
        logger.warning(f"User creation failed: email {email} already exists")
        return None
    except Exception as e:
        logger.error(f"Failed to create user {email}: {e}")
        return None


def get_user_by_email(email: str) -> Optional[Dict[str, Any]]:
    """
    Get user by email address.

    Args:
        email: Email address to look up

    Returns:
        User dictionary or None if not found
    """
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("""
                SELECT id, email, password_hash, full_name, oauth_provider, oauth_provider_id,
                       is_verified, subscription_tier, subscription_status, created_at, last_login, is_active
                FROM users
                WHERE email = ? AND is_active = 1
            """, (email.lower(),))

            row = cursor.fetchone()
            if not row:
                return None

            return {
                'id': row['id'],
                'email': row['email'],
                'password_hash': row['password_hash'],
                'full_name': row['full_name'],
                'oauth_provider': row['oauth_provider'],
                'oauth_provider_id': row['oauth_provider_id'],
                'is_verified': bool(row['is_verified']),
                'subscription_tier': row['subscription_tier'],
                'subscription_status': row['subscription_status'],
                'created_at': row['created_at'],
                'last_login': row['last_login'],
                'is_active': bool(row['is_active'])
            }

    except Exception as e:
        logger.error(f"Failed to get user by email {email}: {e}")
        return None


def get_user_by_id(user_id: str) -> Optional[Dict[str, Any]]:
    """
    Get user by ID.

    Args:
        user_id: User ID to look up

    Returns:
        User dictionary or None if not found
    """
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("""
                SELECT id, email, password_hash, full_name, oauth_provider, oauth_provider_id,
                       is_verified, subscription_tier, subscription_status, created_at, last_login, is_active
                FROM users
                WHERE id = ? AND is_active = 1
            """, (user_id,))

            row = cursor.fetchone()
            if not row:
                return None

            return {
                'id': row['id'],
                'email': row['email'],
                'password_hash': row['password_hash'],
                'full_name': row['full_name'],
                'oauth_provider': row['oauth_provider'],
                'oauth_provider_id': row['oauth_provider_id'],
                'is_verified': bool(row['is_verified']),
                'subscription_tier': row['subscription_tier'],
                'subscription_status': row['subscription_status'],
                'created_at': row['created_at'],
                'last_login': row['last_login'],
                'is_active': bool(row['is_active'])
            }

    except Exception as e:
        logger.error(f"Failed to get user by ID {user_id}: {e}")
        return None


def update_user_login(user_id: str) -> bool:
    """
    Update user's last login timestamp.

    Args:
        user_id: User ID

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET last_login = CURRENT_TIMESTAMP
                WHERE id = ?
            """, (user_id,))

        return True

    except Exception as e:
        logger.error(f"Failed to update last login for user {user_id}: {e}")
        return False


def verify_user_email(verification_token: str) -> bool:
    """
    Verify a user's email address using verification token.

    Args:
        verification_token: Email verification token

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET is_verified = 1, verification_token = NULL
                WHERE verification_token = ?
            """, (verification_token,))

            # Check if any rows were updated
            if conn.total_changes > 0:
                logger.info(f"Email verified successfully")
                return True
            else:
                logger.warning(f"Invalid verification token")
                return False

    except Exception as e:
        logger.error(f"Failed to verify email: {e}")
        return False


def update_user_password(user_id: str, new_password_hash: str) -> bool:
    """
    Update user's password.

    Args:
        user_id: User ID
        new_password_hash: New bcrypt hashed password

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET password_hash = ?
                WHERE id = ?
            """, (new_password_hash, user_id))

        logger.info(f"Password updated for user {user_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to update password for user {user_id}: {e}")
        return False


def save_password_reset_token(email: str, reset_token: str, expires_at: datetime) -> bool:
    """
    Save password reset token for a user.

    Args:
        email: User's email address
        reset_token: Password reset token
        expires_at: Token expiration datetime

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET reset_token = ?, reset_token_expires = ?
                WHERE email = ?
            """, (reset_token, expires_at, email.lower()))

        return True

    except Exception as e:
        logger.error(f"Failed to save reset token for {email}: {e}")
        return False


def get_user_by_reset_token(reset_token: str) -> Optional[Dict[str, Any]]:
    """
    Get user by password reset token if not expired.

    Args:
        reset_token: Password reset token

    Returns:
        User dictionary or None if not found/expired
    """
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("""
                SELECT id, email, password_hash, full_name, reset_token_expires
                FROM users
                WHERE reset_token = ? AND is_active = 1
            """, (reset_token,))

            row = cursor.fetchone()
            if not row:
                return None

            # Check if token is expired
            expires_at = datetime.fromisoformat(row['reset_token_expires']) if row['reset_token_expires'] else None
            if expires_at and expires_at < datetime.now():
                logger.warning(f"Reset token expired for user {row['email']}")
                return None

            return {
                'id': row['id'],
                'email': row['email'],
                'password_hash': row['password_hash'],
                'full_name': row['full_name']
            }

    except Exception as e:
        logger.error(f"Failed to get user by reset token: {e}")
        return None


def clear_password_reset_token(user_id: str) -> bool:
    """
    Clear password reset token after successful reset.

    Args:
        user_id: User ID

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET reset_token = NULL, reset_token_expires = NULL
                WHERE id = ?
            """, (user_id,))

        return True

    except Exception as e:
        logger.error(f"Failed to clear reset token for user {user_id}: {e}")
        return False


def update_user_subscription(user_id: str, tier: str, status: str, stripe_customer_id: Optional[str] = None, stripe_subscription_id: Optional[str] = None) -> bool:
    """
    Update user's subscription information.

    Args:
        user_id: User ID
        tier: Subscription tier ('free', 'pro', 'team')
        status: Subscription status ('active', 'cancelled', 'expired')
        stripe_customer_id: Stripe customer ID (optional)
        stripe_subscription_id: Stripe subscription ID (optional)

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET subscription_tier = ?,
                    subscription_status = ?,
                    stripe_customer_id = ?,
                    stripe_subscription_id = ?
                WHERE id = ?
            """, (tier, status, stripe_customer_id, stripe_subscription_id, user_id))

        logger.info(f"Updated subscription for user {user_id}: {tier} - {status}")
        return True

    except Exception as e:
        logger.error(f"Failed to update subscription for user {user_id}: {e}")
        return False


def get_all_users(limit: int = 100, offset: int = 0, search_query: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Get all users with pagination and optional search.

    Args:
        limit: Maximum number of users to return
        offset: Number of users to skip
        search_query: Optional search string for email or name

    Returns:
        List of user dictionaries
    """
    try:
        with get_db_connection() as conn:
            if search_query:
                # Search by email or name
                search_pattern = f"%{search_query}%"
                cursor = conn.execute("""
                    SELECT id, email, full_name, is_verified, subscription_tier,
                           subscription_status, created_at, last_login, is_active
                    FROM users
                    WHERE (email LIKE ? OR full_name LIKE ?) AND is_active = 1
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (search_pattern, search_pattern, limit, offset))
            else:
                cursor = conn.execute("""
                    SELECT id, email, full_name, is_verified, subscription_tier,
                           subscription_status, created_at, last_login, is_active
                    FROM users
                    WHERE is_active = 1
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (limit, offset))

            users = []
            for row in cursor.fetchall():
                users.append({
                    'id': row['id'],
                    'email': row['email'],
                    'full_name': row['full_name'],
                    'is_verified': bool(row['is_verified']),
                    'subscription_tier': row['subscription_tier'],
                    'subscription_status': row['subscription_status'],
                    'created_at': row['created_at'],
                    'last_login': row['last_login'],
                    'is_active': bool(row['is_active'])
                })

            return users

    except Exception as e:
        logger.error(f"Failed to get users: {e}")
        return []


def get_admin_statistics() -> Dict[str, Any]:
    """
    Get comprehensive admin statistics.

    Returns:
        Dictionary with user statistics
    """
    try:
        with get_db_connection() as conn:
            # Total users
            cursor = conn.execute("SELECT COUNT(*) FROM users WHERE is_active = 1")
            total_users = cursor.fetchone()[0]

            # Verified users
            cursor = conn.execute("SELECT COUNT(*) FROM users WHERE is_verified = 1 AND is_active = 1")
            verified_users = cursor.fetchone()[0]

            # Unverified users
            cursor = conn.execute("SELECT COUNT(*) FROM users WHERE is_verified = 0 AND is_active = 1")
            unverified_users = cursor.fetchone()[0]

            # New users in last 24 hours
            cursor = conn.execute("""
                SELECT COUNT(*) FROM users
                WHERE created_at >= datetime('now', '-1 day') AND is_active = 1
            """)
            new_users_24h = cursor.fetchone()[0]

            # New users in last 7 days
            cursor = conn.execute("""
                SELECT COUNT(*) FROM users
                WHERE created_at >= datetime('now', '-7 days') AND is_active = 1
            """)
            new_users_7d = cursor.fetchone()[0]

            # New users in last 30 days
            cursor = conn.execute("""
                SELECT COUNT(*) FROM users
                WHERE created_at >= datetime('now', '-30 days') AND is_active = 1
            """)
            new_users_30d = cursor.fetchone()[0]

            # Users by subscription tier
            cursor = conn.execute("""
                SELECT subscription_tier, COUNT(*) as count
                FROM users
                WHERE is_active = 1
                GROUP BY subscription_tier
            """)
            users_by_tier = {row['subscription_tier']: row['count'] for row in cursor.fetchall()}

            # Total conversations
            cursor = conn.execute("SELECT COUNT(*) FROM conversations WHERE is_active = 1")
            total_conversations = cursor.fetchone()[0]

            # Total messages
            cursor = conn.execute("SELECT COUNT(*) FROM messages")
            total_messages = cursor.fetchone()[0]

            # Recent registrations (last 10)
            cursor = conn.execute("""
                SELECT email, full_name, created_at, is_verified
                FROM users
                WHERE is_active = 1
                ORDER BY created_at DESC
                LIMIT 10
            """)
            recent_registrations = []
            for row in cursor.fetchall():
                recent_registrations.append({
                    'email': row['email'],
                    'full_name': row['full_name'],
                    'created_at': row['created_at'],
                    'is_verified': bool(row['is_verified'])
                })

            return {
                'total_users': total_users,
                'verified_users': verified_users,
                'unverified_users': unverified_users,
                'new_users_24h': new_users_24h,
                'new_users_7d': new_users_7d,
                'new_users_30d': new_users_30d,
                'users_by_tier': users_by_tier,
                'total_conversations': total_conversations,
                'total_messages': total_messages,
                'recent_registrations': recent_registrations
            }

    except Exception as e:
        logger.error(f"Failed to get admin statistics: {e}")
        return {
            'total_users': 0,
            'verified_users': 0,
            'unverified_users': 0,
            'new_users_24h': 0,
            'new_users_7d': 0,
            'new_users_30d': 0,
            'users_by_tier': {},
            'total_conversations': 0,
            'total_messages': 0,
            'recent_registrations': []
        }


def manually_verify_user(user_id: str) -> bool:
    """
    Manually verify a user (admin action).

    Args:
        user_id: User ID to verify

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET is_verified = 1,
                    verification_token = NULL
                WHERE id = ?
            """, (user_id,))

        logger.info(f"Manually verified user {user_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to manually verify user {user_id}: {e}")
        return False


def create_oauth_user(
    email: str,
    full_name: Optional[str],
    oauth_provider: str,
    oauth_provider_id: str
) -> Optional[str]:
    """
    Create a new user account via OAuth.

    Args:
        email: User's email address from OAuth provider
        full_name: User's full name from OAuth provider
        oauth_provider: OAuth provider ('google' or 'apple')
        oauth_provider_id: Unique user ID from OAuth provider

    Returns:
        User ID if successful, None otherwise
    """
    try:
        user_id = str(uuid.uuid4())
        with get_db_connection() as conn:
            conn.execute("""
                INSERT INTO users (id, email, password_hash, full_name, oauth_provider,
                                   oauth_provider_id, is_verified, verification_token)
                VALUES (?, ?, NULL, ?, ?, ?, 1, NULL)
            """, (user_id, email.lower(), full_name, oauth_provider, oauth_provider_id))

        logger.info(f"Created OAuth user via {oauth_provider}: {email}")
        return user_id

    except sqlite3.IntegrityError as e:
        logger.warning(f"OAuth user creation failed: {e}")
        return None
    except Exception as e:
        logger.error(f"Failed to create OAuth user {email}: {e}")
        return None


def get_user_by_oauth(oauth_provider: str, oauth_provider_id: str) -> Optional[Dict[str, Any]]:
    """
    Get user by OAuth provider and provider ID.

    Args:
        oauth_provider: OAuth provider ('google' or 'apple')
        oauth_provider_id: Unique user ID from OAuth provider

    Returns:
        User dictionary or None if not found
    """
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("""
                SELECT id, email, password_hash, full_name, oauth_provider, oauth_provider_id,
                       is_verified, subscription_tier, subscription_status, created_at,
                       last_login, is_active
                FROM users
                WHERE oauth_provider = ? AND oauth_provider_id = ? AND is_active = 1
            """, (oauth_provider, oauth_provider_id))

            row = cursor.fetchone()
            if not row:
                return None

            return {
                'id': row['id'],
                'email': row['email'],
                'password_hash': row['password_hash'],
                'full_name': row['full_name'],
                'oauth_provider': row['oauth_provider'],
                'oauth_provider_id': row['oauth_provider_id'],
                'is_verified': bool(row['is_verified']),
                'subscription_tier': row['subscription_tier'],
                'subscription_status': row['subscription_status'],
                'created_at': row['created_at'],
                'last_login': row['last_login'],
                'is_active': bool(row['is_active'])
            }

    except Exception as e:
        logger.error(f"Failed to get user by OAuth {oauth_provider}/{oauth_provider_id}: {e}")
        return None


def link_oauth_to_existing_user(
    user_id: str,
    oauth_provider: str,
    oauth_provider_id: str
) -> bool:
    """
    Link OAuth credentials to an existing user account.

    Args:
        user_id: Existing user ID
        oauth_provider: OAuth provider ('google' or 'apple')
        oauth_provider_id: Unique user ID from OAuth provider

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("""
                UPDATE users
                SET oauth_provider = ?,
                    oauth_provider_id = ?,
                    is_verified = 1
                WHERE id = ?
            """, (oauth_provider, oauth_provider_id, user_id))

        logger.info(f"Linked {oauth_provider} OAuth to user {user_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to link OAuth to user {user_id}: {e}")
        return False


def migrate_database_for_oauth() -> bool:
    """
    Migrate existing database to support OAuth.
    Safe to run multiple times - only adds columns if they don't exist.

    Returns:
        True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            # Check if oauth columns exist
            cursor = conn.execute("PRAGMA table_info(users)")
            columns = [row['name'] for row in cursor.fetchall()]

            if 'oauth_provider' not in columns:
                logger.info("Adding oauth_provider column to users table")
                conn.execute("""
                    ALTER TABLE users
                    ADD COLUMN oauth_provider TEXT CHECK(oauth_provider IN ('google', 'apple', NULL))
                """)

            if 'oauth_provider_id' not in columns:
                logger.info("Adding oauth_provider_id column to users table")
                conn.execute("""
                    ALTER TABLE users
                    ADD COLUMN oauth_provider_id TEXT
                """)

            # Check if password_hash has NOT NULL constraint (would prevent OAuth users)
            if 'password_hash' in columns:
                cursor = conn.execute("""
                    SELECT sql FROM sqlite_master
                    WHERE type='table' AND name='users'
                """)
                table_sql = cursor.fetchone()
                if table_sql and 'password_hash TEXT NOT NULL' in table_sql['sql']:
                    # Note: SQLite doesn't support ALTER COLUMN to drop NOT NULL
                    # OAuth users will need password_hash to be nullable
                    # If this is a fresh install, the schema already has password_hash as nullable
                    logger.info("Note: Existing users table has NOT NULL constraint on password_hash")
                    logger.info("OAuth sign-in requires nullable password_hash field")
                    logger.info("For fresh installs, this is already handled by the current schema")

        logger.info("OAuth migration completed successfully")
        return True

    except Exception as e:
        logger.error(f"Failed to migrate database for OAuth: {e}")
        return False


# ====== Bookmark Functions ======

def save_bookmark(
    user_session_id: str,
    query: str,
    answer: str,
    references: Optional[List[Dict]] = None,
    num_papers: int = 0,
    evidence_strength: Optional[str] = None,
    user_id: Optional[str] = None
) -> Optional[str]:
    """
    Save a bookmark to the database.

    Args:
        user_session_id: Session ID (for guest users)
        query: The original query
        answer: The assistant's answer
        references: List of paper references
        num_papers: Number of papers cited
        evidence_strength: Evidence quality ('high', 'moderate', 'low')
        user_id: User ID (for authenticated users)

    Returns:
        str: Bookmark ID if successful, None otherwise
    """
    try:
        bookmark_id = str(uuid.uuid4())
        references_json = json.dumps(references) if references else None

        with get_db_connection() as conn:
            conn.execute("""
                INSERT INTO bookmarks
                (id, user_session_id, user_id, query, answer, references, num_papers, evidence_strength)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (bookmark_id, user_session_id, user_id, query, answer, references_json, num_papers, evidence_strength))

        logger.info(f"Saved bookmark {bookmark_id} for user/session {user_id or user_session_id}")
        return bookmark_id

    except Exception as e:
        logger.error(f"Failed to save bookmark: {e}")
        return None


def get_bookmarks(
    user_session_id: Optional[str] = None,
    user_id: Optional[str] = None,
    limit: int = 50,
    offset: int = 0
) -> List[Dict]:
    """
    Retrieve bookmarks for a user or session.

    Args:
        user_session_id: Session ID (for guest users)
        user_id: User ID (for authenticated users)
        limit: Maximum number of bookmarks to return
        offset: Number of bookmarks to skip

    Returns:
        List of bookmark dictionaries
    """
    try:
        with get_db_connection() as conn:
            cursor = conn.cursor()

            if user_id:
                cursor.execute("""
                    SELECT id, query, answer, paper_references, num_papers, evidence_strength, created_at
                    FROM bookmarks
                    WHERE user_id = ?
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (user_id, limit, offset))
            elif user_session_id:
                cursor.execute("""
                    SELECT id, query, answer, paper_references, num_papers, evidence_strength, created_at
                    FROM bookmarks
                    WHERE user_session_id = ?
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (user_session_id, limit, offset))
            else:
                return []

            rows = cursor.fetchall()
            bookmarks = []

            for row in rows:
                bookmark = {
                    'id': row['id'],
                    'query': row['query'],
                    'answer': row['answer'],
                    'references': json.loads(row['paper_references']) if row['paper_references'] else [],
                    'num_papers': row['num_papers'],
                    'evidence_strength': row['evidence_strength'],
                    'timestamp': row['created_at']
                }
                bookmarks.append(bookmark)

            return bookmarks

    except Exception as e:
        logger.error(f"Failed to retrieve bookmarks: {e}")
        return []


def delete_bookmark(bookmark_id: str) -> bool:
    """
    Delete a bookmark by ID.

    Args:
        bookmark_id: The bookmark ID to delete

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with get_db_connection() as conn:
            conn.execute("DELETE FROM bookmarks WHERE id = ?", (bookmark_id,))

        logger.info(f"Deleted bookmark {bookmark_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to delete bookmark {bookmark_id}: {e}")
        return False


# ====== Shared Link Functions ======

def create_shared_link(
    query: str,
    answer: str,
    references: Optional[List[Dict]] = None,
    expires_days: int = 30
) -> Optional[str]:
    """
    Create a shareable link for a response.

    Args:
        query: The original query
        answer: The assistant's answer
        references: List of paper references
        expires_days: Number of days until link expires (default: 30)

    Returns:
        str: Share ID if successful, None otherwise
    """
    try:
        link_id = str(uuid.uuid4())
        share_id = str(uuid.uuid4())[:12]  # Short ID for URLs
        references_json = json.dumps(references) if references else None

        from datetime import datetime, timedelta
        expires_at = datetime.now() + timedelta(days=expires_days)

        with get_db_connection() as conn:
            conn.execute("""
                INSERT INTO shared_links
                (id, share_id, query, answer, references, expires_at)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (link_id, share_id, query, answer, references_json, expires_at))

        logger.info(f"Created shared link {share_id} (expires: {expires_at})")
        return share_id

    except Exception as e:
        logger.error(f"Failed to create shared link: {e}")
        return None


def get_shared_link(share_id: str) -> Optional[Dict]:
    """
    Retrieve a shared link by share ID.

    Args:
        share_id: The share ID

    Returns:
        Dict containing the shared data, or None if not found/expired
    """
    try:
        from datetime import datetime

        with get_db_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT share_id, query, answer, paper_references, created_at, expires_at
                FROM shared_links
                WHERE share_id = ?
            """, (share_id,))

            row = cursor.fetchone()

            if not row:
                logger.warning(f"Shared link not found: {share_id}")
                return None

            # Check if expired
            expires_at = datetime.fromisoformat(row['expires_at'])
            if datetime.now() > expires_at:
                logger.info(f"Shared link expired: {share_id}")
                # Delete expired link
                conn.execute("DELETE FROM shared_links WHERE share_id = ?", (share_id,))
                return None

            shared_data = {
                'share_id': row['share_id'],
                'query': row['query'],
                'answer': row['answer'],
                'references': json.loads(row['paper_references']) if row['paper_references'] else [],
                'timestamp': row['created_at'],
                'expires': row['expires_at']
            }

            return shared_data

    except Exception as e:
        logger.error(f"Failed to retrieve shared link {share_id}: {e}")
        return None


def cleanup_expired_shared_links() -> int:
    """
    Delete all expired shared links.

    Returns:
        int: Number of links deleted
    """
    try:
        from datetime import datetime

        with get_db_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                DELETE FROM shared_links
                WHERE expires_at < ?
            """, (datetime.now(),))

            deleted_count = cursor.rowcount

        if deleted_count > 0:
            logger.info(f"Cleaned up {deleted_count} expired shared links")

        return deleted_count

    except Exception as e:
        logger.error(f"Failed to cleanup expired shared links: {e}")
        return 0


if __name__ == '__main__':
    """
    Test the database module independently.
    Run with: python database.py
    """
    print("Testing database module...")
    print(f"Database path: {DB_PATH}")

    # Initialize database
    print("\n1. Initializing database...")
    if init_db():
        print("✓ Database initialized successfully")
    else:
        print("✗ Database initialization failed")
        exit(1)

    # Create test conversation
    print("\n2. Creating test conversation...")
    test_conv_id = str(uuid.uuid4())
    test_user_id = "test_user_123"
    if save_conversation(test_conv_id, test_user_id, "Test conversation"):
        print(f"✓ Created conversation: {test_conv_id}")
    else:
        print("✗ Failed to create conversation")
        exit(1)

    # Save test messages
    print("\n3. Saving test messages...")
    if save_message(test_conv_id, "user", "What is propofol?"):
        print("✓ Saved user message")
    else:
        print("✗ Failed to save user message")
        exit(1)

    if save_message(
        test_conv_id,
        "assistant",
        "Propofol is an intravenous anesthetic...",
        references=[{"title": "Test paper", "pmid": "12345"}],
        num_papers=1,
        evidence_strength="high"
    ):
        print("✓ Saved assistant message")
    else:
        print("✗ Failed to save assistant message")
        exit(1)

    # Retrieve conversation
    print("\n4. Retrieving conversation...")
    conv = get_conversation(test_conv_id)
    if conv:
        print(f"✓ Retrieved conversation: {conv['title']}")
        print(f"  Messages: {len(conv['messages'])}")
    else:
        print("✗ Failed to retrieve conversation")
        exit(1)

    # Get conversations list
    print("\n5. Getting conversations list...")
    convs = get_conversations(test_user_id)
    if convs:
        print(f"✓ Retrieved {len(convs)} conversations")
    else:
        print("✗ Failed to get conversations list")
        exit(1)

    # Get database stats
    print("\n6. Getting database stats...")
    stats = get_database_stats()
    print(f"✓ Database stats:")
    print(f"  Path: {stats['database_path']}")
    print(f"  Size: {stats.get('database_size_mb', 0)} MB")
    print(f"  Conversations: {stats.get('total_conversations', 0)}")
    print(f"  Messages: {stats.get('total_messages', 0)}")

    # Clean up test data
    print("\n7. Cleaning up test data...")
    if delete_conversation(test_conv_id, hard_delete=True):
        print("✓ Deleted test conversation")
    else:
        print("✗ Failed to delete test conversation")

    print("\n✅ All tests passed!")

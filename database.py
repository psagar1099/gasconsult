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
    logger.warning(f"Production DB directory {DB_DIR} doesn't exist, using fallback")
    DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)
    logger.info(f"Using fallback database path: {DB_PATH}")


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
            # Create conversations table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS conversations (
                    id TEXT PRIMARY KEY,
                    user_session_id TEXT NOT NULL,
                    title TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    message_count INTEGER DEFAULT 0,
                    is_active BOOLEAN DEFAULT 1
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
                CREATE INDEX IF NOT EXISTS idx_conversations_session
                ON conversations(user_session_id, created_at DESC)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_conversations_active
                ON conversations(user_session_id, is_active, created_at DESC)
            """)

            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_messages_conversation
                ON messages(conversation_id, created_at ASC)
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

        with get_db_connection() as conn:
            # Insert message
            conn.execute("""
                INSERT INTO messages
                (id, conversation_id, role, content, paper_references, num_papers, evidence_strength)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (message_id, conversation_id, role, content, references_json, num_papers, evidence_strength))

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
                        message['evidence_strength'] = msg_row['evidence_strength']

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

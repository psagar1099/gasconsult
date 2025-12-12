#!/usr/bin/env python3
"""
Database migration script to add OAuth columns to existing users table.
This fixes the "no such column: oauth_provider" error.

Run this script once to update the production database.
"""

import sys
import os

# Add current directory to path to import database module
sys.path.insert(0, os.path.dirname(__file__))

from database import migrate_database_for_oauth, get_db_connection, DB_PATH, logger
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    print("=" * 60)
    print("OAuth Database Migration Script")
    print("=" * 60)
    print(f"\nDatabase path: {DB_PATH}")
    print(f"Database exists: {os.path.exists(DB_PATH)}")

    if not os.path.exists(DB_PATH):
        print("\n❌ ERROR: Database does not exist!")
        print("Please ensure the database has been initialized first.")
        return 1

    print("\n📊 Checking current table structure...")
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("PRAGMA table_info(users)")
            columns = [row['name'] for row in cursor.fetchall()]

            print(f"\nCurrent columns in 'users' table:")
            for col in columns:
                print(f"  ✓ {col}")

            # Check if migration is needed
            needs_migration = False
            if 'oauth_provider' not in columns:
                print("\n⚠️  Missing column: oauth_provider")
                needs_migration = True
            if 'oauth_provider_id' not in columns:
                print("⚠️  Missing column: oauth_provider_id")
                needs_migration = True

            if not needs_migration:
                print("\n✅ Database already has OAuth columns - no migration needed!")
                return 0

    except Exception as e:
        print(f"\n❌ Error checking database: {e}")
        return 1

    print("\n🔧 Running migration...")
    if migrate_database_for_oauth():
        print("\n✅ Migration completed successfully!")

        # Verify the migration
        print("\n📊 Verifying migration...")
        try:
            with get_db_connection() as conn:
                cursor = conn.execute("PRAGMA table_info(users)")
                columns = [row['name'] for row in cursor.fetchall()]

                if 'oauth_provider' in columns and 'oauth_provider_id' in columns:
                    print("✅ OAuth columns successfully added!")
                    print("\nUpdated columns in 'users' table:")
                    for col in columns:
                        if col in ['oauth_provider', 'oauth_provider_id']:
                            print(f"  ✓ {col} (NEW)")
                        else:
                            print(f"  ✓ {col}")

                    # Count users
                    cursor = conn.execute("SELECT COUNT(*) as count FROM users WHERE is_active = 1")
                    user_count = cursor.fetchone()['count']
                    print(f"\n👥 Total active users: {user_count}")

                    return 0
                else:
                    print("❌ Migration verification failed - columns not found")
                    return 1

        except Exception as e:
            print(f"❌ Error verifying migration: {e}")
            return 1
    else:
        print("\n❌ Migration failed!")
        print("Please check the error logs above for details.")
        return 1

if __name__ == '__main__':
    exit_code = main()
    print("\n" + "=" * 60)
    sys.exit(exit_code)

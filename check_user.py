#!/usr/bin/env python3
"""
Diagnostic script to check user account status.
Usage: python check_user.py <email>
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from database import get_user_by_email, get_db_connection, DB_PATH

def check_user(email):
    print("=" * 70)
    print(f"User Account Diagnostic: {email}")
    print("=" * 70)
    print(f"\nDatabase path: {DB_PATH}")
    print(f"Database exists: {os.path.exists(DB_PATH)}")

    if not os.path.exists(DB_PATH):
        print("\n❌ ERROR: Database does not exist!")
        return 1

    print("\n📊 Checking database schema...")
    try:
        with get_db_connection() as conn:
            cursor = conn.execute("PRAGMA table_info(users)")
            columns = [row['name'] for row in cursor.fetchall()]

            print("\nColumns in users table:")
            required_cols = ['id', 'email', 'password_hash', 'oauth_provider', 'oauth_provider_id', 'is_verified']
            for col in columns:
                marker = "✓" if col in required_cols else " "
                print(f"  {marker} {col}")

            # Check for OAuth columns
            has_oauth = 'oauth_provider' in columns and 'oauth_provider_id' in columns
            print(f"\nOAuth columns present: {'✅ Yes' if has_oauth else '❌ No'}")

    except Exception as e:
        print(f"\n❌ Error checking schema: {e}")
        return 1

    print(f"\n👤 Looking up user: {email}")
    print("-" * 70)

    try:
        user_data = get_user_by_email(email)

        if not user_data:
            print(f"\n❌ User not found: {email}")

            # Check if user exists but is inactive
            with get_db_connection() as conn:
                cursor = conn.execute(
                    "SELECT id, email, is_active FROM users WHERE email = ?",
                    (email.lower(),)
                )
                row = cursor.fetchone()
                if row:
                    print(f"   Found inactive user: is_active = {row['is_active']}")
                else:
                    print("   User does not exist in database at all")

            return 1

        print(f"\n✅ User found!")
        print("\nAccount Details:")
        print(f"  ID: {user_data['id']}")
        print(f"  Email: {user_data['email']}")
        print(f"  Full Name: {user_data.get('full_name', 'N/A')}")
        print(f"  Is Verified: {user_data.get('is_verified', False)}")
        print(f"  Is Active: {user_data.get('is_active', False)}")
        print(f"  Subscription Tier: {user_data.get('subscription_tier', 'N/A')}")
        print(f"  Created At: {user_data.get('created_at', 'N/A')}")
        print(f"  Last Login: {user_data.get('last_login', 'N/A')}")

        print("\nAuthentication Details:")
        has_password = bool(user_data.get('password_hash'))
        has_oauth = bool(user_data.get('oauth_provider'))

        print(f"  Password Hash: {'✅ Set' if has_password else '❌ NULL'}")
        print(f"  OAuth Provider: {user_data.get('oauth_provider') or 'None'}")
        print(f"  OAuth Provider ID: {user_data.get('oauth_provider_id') or 'None'}")

        print("\nAccount Type:")
        if has_password and not has_oauth:
            print("  ✅ Email/Password Account (Normal)")
        elif has_oauth and not has_password:
            print("  ⚠️  OAuth-Only Account (Google/Apple)")
            print(f"     → Must use {user_data['oauth_provider'].title()} Sign In button")
        elif has_password and has_oauth:
            print("  🔗 Linked Account (Password + OAuth)")
            print("     → Can login with either method")
        else:
            print("  ❌ BROKEN ACCOUNT (No password AND no OAuth)")
            print("     → This account cannot login!")

        print("\nLogin Status:")
        if not user_data.get('is_active'):
            print("  ❌ Account is DEACTIVATED")
        elif not has_password and not has_oauth:
            print("  ❌ Cannot login (no authentication method)")
        elif has_password:
            print("  ✅ Can login with email/password")
        else:
            print(f"  ✅ Can login with {user_data['oauth_provider'].title()} Sign In")

        return 0

    except Exception as e:
        print(f"\n❌ Error looking up user: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python check_user.py <email>")
        print("\nExample:")
        print("  python check_user.py sagar@gasconsult.ai")
        sys.exit(1)

    email = sys.argv[1]
    exit_code = check_user(email)

    print("\n" + "=" * 70)
    sys.exit(exit_code)

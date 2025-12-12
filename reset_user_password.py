#!/usr/bin/env python3
"""
Manual password reset script for emergency account recovery.
Usage: python reset_user_password.py <email> <new_password>

IMPORTANT: Only use this for account recovery.
Users should normally use the "Forgot Password" feature.
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from database import get_user_by_email, update_user_password, get_db_connection
import bcrypt

def reset_password(email, new_password):
    print("=" * 70)
    print(f"Manual Password Reset: {email}")
    print("=" * 70)

    if len(new_password) < 8:
        print("\n❌ ERROR: Password must be at least 8 characters long")
        return 1

    print(f"\n👤 Looking up user: {email}")
    try:
        user_data = get_user_by_email(email)

        if not user_data:
            print(f"❌ User not found: {email}")
            return 1

        user_id = user_data['id']
        print(f"✅ User found: {user_id}")

        # Check if it's an OAuth-only account
        if user_data.get('oauth_provider') and not user_data.get('password_hash'):
            print(f"\n⚠️  WARNING: This is an OAuth-only account ({user_data['oauth_provider']})")
            print("Setting a password will convert it to a hybrid account.")
            response = input("Continue? (yes/no): ")
            if response.lower() != 'yes':
                print("Cancelled.")
                return 0

        # Hash the new password
        print("\n🔒 Hashing new password...")
        password_hash = bcrypt.hashpw(new_password.encode('utf-8'), bcrypt.gensalt()).decode('utf-8')

        # Update password in database
        print("💾 Updating password in database...")
        if update_user_password(user_id, password_hash):
            print("✅ Password updated successfully!")
            print(f"\nUser {email} can now log in with the new password.")
            return 0
        else:
            print("❌ Failed to update password in database")
            return 1

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python reset_user_password.py <email> <new_password>")
        print("\nExample:")
        print('  python reset_user_password.py sagar@gasconsult.ai "NewSecurePass123!"')
        print("\nSecurity Note:")
        print("  - Password must be at least 8 characters")
        print("  - Use quotes around password if it contains spaces or special characters")
        sys.exit(1)

    email = sys.argv[1]
    new_password = sys.argv[2]

    exit_code = reset_password(email, new_password)

    print("\n" + "=" * 70)
    sys.exit(exit_code)

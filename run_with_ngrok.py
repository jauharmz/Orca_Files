#!/usr/bin/env python3
"""
Run Flask app with ngrok tunnel for Google Colab.
No IP confirmation needed - direct access to UI.
"""

from pyngrok import ngrok
import subprocess
import sys
import time
import threading

def main():
    print("=" * 70)
    print("🧪 ORCA OUTPUT PARSER - NGROK TUNNEL")
    print("=" * 70)
    print()

    # Start ngrok tunnel
    print("🔗 Creating ngrok tunnel...")
    public_url = ngrok.connect(5000)
    print(f"✅ Your app is accessible at: {public_url}")
    print()
    print("=" * 70)
    print("🌐 Starting Flask server on port 5000...")
    print("=" * 70)
    print()

    # Start Flask (this will block)
    subprocess.run([sys.executable, "app.py"])

if __name__ == "__main__":
    main()

"""
Automated Flask + Localtunnel Runner for Google Colab

This script automatically starts Flask and creates a public tunnel.
Displays your IP address for easy copy/paste into localtunnel landing page.

Usage:
    python run_with_tunnel.py
"""

import subprocess
import threading
import time
import sys
import urllib.request

def get_public_ip():
    """Fetch public IP address."""
    try:
        with urllib.request.urlopen('https://api.ipify.org', timeout=5) as response:
            return response.read().decode('utf-8')
    except:
        try:
            with urllib.request.urlopen('https://ifconfig.me', timeout=5) as response:
                return response.read().decode('utf-8')
        except:
            return "Unable to fetch (check manually at https://api.ipify.org)"

def run_flask():
    """Run Flask server."""
    print("\n" + "="*70)
    print("🌐 Starting Flask server on port 5000...")
    print("="*70)
    subprocess.run([sys.executable, "app.py"])

def run_tunnel():
    """Run localtunnel after Flask starts."""
    print("\n" + "="*70)
    print("⏳ Waiting 3 seconds for Flask to start...")
    print("="*70)
    time.sleep(3)

    # Fetch IP address
    print("\n" + "="*70)
    print("🔍 Detecting your public IP address...")
    print("="*70)
    ip_address = get_public_ip()

    print("\n" + "="*70)
    print("📋 YOUR IP ADDRESS (copy this!):")
    print("="*70)
    print(f"\n    {ip_address}\n")
    print("="*70)
    print("⚠️  IMPORTANT: When you open the tunnel URL, you'll see a landing page.")
    print("   Just paste this IP address and click Continue!")
    print("="*70)

    print("\n" + "="*70)
    print("🔗 Creating public tunnel with localtunnel...")
    print("="*70)
    print("Your app will be accessible at: https://XXXXX.loca.lt")
    print("(The URL will appear below)")
    print("="*70 + "\n")

    subprocess.run(["lt", "--port", "5000"])

def main():
    print("\n" + "="*70)
    print("🚀 ORCA PARSER - AUTOMATED LAUNCH WITH TUNNEL")
    print("="*70)
    print("\nThis will:")
    print("  1. Start Flask server on port 5000")
    print("  2. Create public tunnel via localtunnel")
    print("  3. Give you a public URL to access the UI")
    print("\n" + "="*70)
    
    # Start Flask in a thread
    flask_thread = threading.Thread(target=run_flask, daemon=True)
    flask_thread.start()
    
    # Run tunnel (this blocks)
    run_tunnel()

if __name__ == "__main__":
    main()

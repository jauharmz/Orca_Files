"""
Automated Flask + Localtunnel Runner for Google Colab

This script automatically starts Flask and creates a public tunnel.

Usage:
    python run_with_tunnel.py
"""

import subprocess
import threading
import time
import sys

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

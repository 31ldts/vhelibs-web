#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VHELIBS Web – entry point.

Usage:
    python run.py [--host HOST] [--port PORT] [--cache-dir PATH] [--no-browser] [--debug]
"""
import argparse
import os
import sys
import webbrowser
import threading


def main():
    p = argparse.ArgumentParser(
        description="VHELIBS Web – Validation Helper for LIgands and Binding Sites"
    )
    p.add_argument("--host",      default="127.0.0.1",
                   help="Bind address (default: 127.0.0.1)")
    p.add_argument("--port",      default=8000, type=int,
                   help="Port number (default: 8000)")
    p.add_argument("--cache-dir", default=None,
                   help="Cache directory for downloaded PDB/EDS files "
                        "(default: ~/.cache/vhelibs)")
    p.add_argument("--no-browser", action="store_true",
                   help="Do not open browser automatically")
    p.add_argument("--debug", action="store_true",
                   help="Enable Flask debug mode (do NOT use in production)")
    args = p.parse_args()

    cache_dir = args.cache_dir or os.path.join(
        os.path.expanduser("~"), ".cache", "vhelibs"
    )
    os.makedirs(cache_dir, exist_ok=True)

    from app import create_app
    app = create_app(cache_dir=cache_dir)

    url = f"http://{args.host}:{args.port}"
    if not args.no_browser and not args.debug:
        # Open browser after a short delay so Flask has time to start
        threading.Timer(1.2, lambda: webbrowser.open(url)).start()

    print(f"VHELIBS Web server starting at {url}")
    print(f"Cache directory: {cache_dir}")
    print("Press Ctrl+C to stop.\n")

    app.run(
        host=args.host,
        port=args.port,
        debug=args.debug,
        use_reloader=False,   # reloader conflicts with the browser-open timer
        threaded=True,        # allow concurrent analysis threads
    )


if __name__ == "__main__":
    main()

"""Static server for the benchmark explorer.

Adds `Cache-Control: no-cache` so browsers revalidate every file (index.html and
data/*.json) instead of masking updates with a stale heuristic cache. Unchanged
files still return a fast 304 via If-Modified-Since. Binds 0.0.0.0 for LAN access.
"""
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import sys

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8848


class Handler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()


if __name__ == "__main__":
    print(f"serving on 0.0.0.0:{PORT}")
    ThreadingHTTPServer(("0.0.0.0", PORT), Handler).serve_forever()

import os
import sys

def debug(*args, **kwargs):
    if not os.getenv("DEBUG", None):
        return
    print(*args, **kwargs, file=sys.stderr)


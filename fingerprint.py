import sys, hashlib
name, h, n = None, None, 0
def out():
    if name: print(f"{name}\t{n}\t{h.hexdigest()}")
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith(">"):
            out(); name, h, n = line[1:].split()[0], hashlib.md5(), 0
        else:
            s = line.strip().upper().encode(); h.update(s); n += len(s)
out()

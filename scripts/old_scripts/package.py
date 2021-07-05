import tarfile
from glob import glob
tar = tarfile.open("../package.tar.gz", "w:gz")
for name in glob("../makefile"):
    tar.add(name)
for name in glob("../*.c"):
    tar.add(name)
for name in glob("../*.h"):
    tar.add(name)
for name in glob("../rng/*.c"):
    tar.add(name)
for name in glob("../rng/*.h"):
    tar.add(name)
for name in glob("./*.sh"):
    tar.add(name)
tar.close()

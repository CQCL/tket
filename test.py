import os
import subprocess
import re

f = open("blah.txt", "w")
subprocess.call(["cat", "testsssssss"], stdout=f)
f.close()
with open("blah.txt", "r") as hef:
    shimgen_help = hef.read()
    ccache_bin = re.search(r"C:\\+ProgramData\\+Chocolatey\\+lib\\+ccache\\+tools\\+ccache-4\.8\.3-windows-x86_64\\+ccache\.exe", shimgen_help).group()
    print(f"Found ccache binary at {ccache_bin}")
#with open(os.environ["GITHUB_OUTPUT"], "a") as f:
#    print("{0}={1}".format("ccache_binary", ccache_bin), file=f)

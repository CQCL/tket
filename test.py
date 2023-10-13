import os
import subprocess
import re

try:
    shimgen_help = subprocess.check_output("cat testsssssss", shell=True, text=True)
except Exception:
    pass
ccache_bin = re.search(r"C:\\+ProgramData\\+Chocolatey\\+lib\\+ccache\\+tools\\+ccache-4\.8\.3-windows-x86_64\\+ccache\.exe", shimgen_help).group()
print(f"Found ccache binary at {ccache_bin}")
#with open(os.environ["GITHUB_OUTPUT"], "a") as f:
#    print("{0}={1}".format("ccache_binary", ccache_bin), file=f)

import subprocess
import re

arch = str(subprocess.check_output("cat testsssssss", shell=True))
print(arch)
arch = re.search(r"C:\\+ProgramData\\+Chocolatey\\+lib\\+ccache\\+tools\\+ccache-4\.8\.3-windows-x86_64\\+ccache\.exe", arch)
print(arch.group())

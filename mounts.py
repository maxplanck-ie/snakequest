#!/usr/bin/env python
import os
import os.path
import subprocess
import shutil
import glob

# Handle automounter updates
if os.path.exists("/export/automount"):
    try:
        os.mkdir("/etc/automount")
    except:
        pass
    for f in glob.glob("/export/automount/*"):
        if os.path.basename(f) in ["auto.master", "yp.conf", "nsswitch.conf"]:
            shutil.copy(f, "/etc/{}".format(os.path.basename(f)))
        else:
            shutil.copy(f, "/etc/automount/{}".format(os.path.basename(f)))
    subprocess.check_call(["/usr/sbin/rsyslogd"])
    #subprocess.check_call(["domainname", "solsys1.immunbio.mpg.de"])
    subprocess.check_call(["/sbin/rpcbind"])
    subprocess.check_call(["/sbin/rpc.statd", "--no-notify"])
    subprocess.check_call(["/usr/sbin/ypbind"])
    subprocess.check_call(["/usr/sbin/automount"])
    subprocess.check_call(["sleep", "15"])


import json, re
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

linLimit = 100.0;
angLimit = 5000.0;

params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
plt.rcParams.update(params)

# Remove comments in json string before parsing
config = json.loads(re.sub("//.*","", \
        open("input_case9218Jietal2015.json").read(),flags=re.MULTILINE))
settings = config["simulation"]
tMax =settings["maximum-time"]*1000.0

xt = np.array(settings["linear-acceleration"]["xt"])
xv = np.array(settings["linear-acceleration"]["xv"])
yt = np.array(settings["linear-acceleration"]["yt"])
yv = np.array(settings["linear-acceleration"]["yv"])
zt = np.array(settings["linear-acceleration"]["zt"])
zv = np.array(settings["linear-acceleration"]["zv"])

plt.plot(xt, xv, 'r', linewidth=4.0, label='X');
plt.plot(yt, yv, 'b', linewidth=4.0, label='Y');
plt.plot(zt, zv, 'k', linewidth=4.0, label='Z');
plt.xlabel('Time (ms)')
plt.ylabel('Linear Acceleration (g)')
plt.legend(loc='best')

plt.xlim(0, tMax)
plt.ylim(-linLimit, linLimit)
plt.savefig('linearAcceleration.png')

xt = np.array(settings["angular-acceleration"]["xt"])
xv = np.array(settings["angular-acceleration"]["xv"])
yt = np.array(settings["angular-acceleration"]["yt"])
yv = np.array(settings["angular-acceleration"]["yv"])
zt = np.array(settings["angular-acceleration"]["zt"])
zv = np.array(settings["angular-acceleration"]["zv"])

plt.figure()
plt.plot(xt, xv, 'r', linewidth=4.0, label='X');
plt.plot(yt, yv, 'b', linewidth=4.0, label='Y');
plt.plot(zt, zv, 'k', linewidth=4.0, label='Z');
plt.xlabel('Time (ms)')
plt.ylabel(r'Angular Acceleration (rad/$s^2$)')
plt.legend(loc='best')

plt.xlim(0, tMax)
plt.ylim(-angLimit, angLimit)
plt.savefig('angularAcceleration.png')

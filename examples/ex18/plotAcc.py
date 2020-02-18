import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Validation-case from FemTech mapping
# x <--> -z
# y <--> +x
# z <--> -y

linLimit = 100.0;
angLimit = 5000.0;

params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
plt.rcParams.update(params)

config = json.loads(open("input_case9218Jietal2015.json").read())
settings = config["simulation"]
tMax =settings["maximum-time"]*1000.0

xtF = np.array(settings["linear-acceleration"]["xt"])
xvF = np.array(settings["linear-acceleration"]["xv"])
ytF = np.array(settings["linear-acceleration"]["yt"])
yvF = np.array(settings["linear-acceleration"]["yv"])
ztF = np.array(settings["linear-acceleration"]["zt"])
zvF = np.array(settings["linear-acceleration"]["zv"])

#Map to validation case values
xt = ztF; xv = -zvF;
yt = xtF; yv = xvF;
zt = ytF; zv = -yvF;

plt.plot(xt, xv, 'r', linewidth=4.0, label='X');
plt.plot(yt, yv, 'b', linewidth=4.0, label='Y');
plt.plot(zt, zv, 'k', linewidth=4.0, label='Z');
plt.xlabel('Time (ms)')
plt.ylabel('Linear Acceleration (g)')
plt.legend(loc='best')

plt.xlim(0, tMax)
plt.ylim(-linLimit, linLimit)
plt.savefig('linearAcceleration.png')

xtF = np.array(settings["angular-acceleration"]["xt"])
xvF = np.array(settings["angular-acceleration"]["xv"])
ytF = np.array(settings["angular-acceleration"]["yt"])
yvF = np.array(settings["angular-acceleration"]["yv"])
ztF = np.array(settings["angular-acceleration"]["zt"])
zvF = np.array(settings["angular-acceleration"]["zv"])

#Map to validation case values
xt = ztF; xv = -zvF;
yt = xtF; yv = xvF;
zt = ytF; zv = -yvF;

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

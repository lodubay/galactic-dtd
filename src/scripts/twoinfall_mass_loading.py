import matplotlib.pyplot as plt
import vice
import paths

output = '../data/migration/diffusion/twoinfall/exponential_timescale15.vice'
zones = [4, 6, 8, 10, 12, 14]

fig, ax = plt.subplots()

for z in zones:
    hist = vice.history(output + '/zone%d' % (z * 10))
    ax.plot(hist['time'][1:], hist['ofr'][1:], label=hist['eta_0'][1])

ax.set_xlabel('Time [Gyr]')
ax.set_ylabel('OFR')
ax.legend()
plt.savefig(paths.figures / 'twoinfall_mass_loading.png', dpi=300)
import numpy as np, matplotlib.pyplot as plt
from astropy import units as u, constants as c
from astropy.coordinates import SkyCoord, AltAz, concatenate, EarthLocation
from astropy.time import Time

# open data file
#data = np.load('data_3C273_fixed.npz')
data = np.load('Downloads/DSA_lab/data_3C273_fixed.npz')

I_0 = 42*u.Jy # from NED, at 1.41e9 Hz

dec_deg = data['DEC_degrees']*u.deg
nu = data['frequencies_Hz']*u.Hz

# calculate s

ovro = EarthLocation(lat=37.2317*u.deg,lon=-118.2951*u.deg)
time = Time(data['times_MJD'], format = 'mjd')
aa = AltAz(obstime = time, location=ovro)

qso = SkyCoord(ra = 187.277915*u.deg, dec = 2.052388*u.deg)

qso_aa = qso.transform_to(aa)

s = [qso_aa.az, qso_aa.alt]

s_proj = [np.sin(s[0]), np.cos(s[0]), np.sin(s[1])]

s_unit_proj = s_proj/np.linalg.norm(s_proj)


def g_ab(baseline):

    base = baseline
    # time = 0

    u_m = list(data['u_meters'][:,base]) #ntimes, nbaselines
    v_m = list(data['v_meters'][:,base]) #ntimes, nbaselines

    b = np.array([u_m, v_m, list(np.zeros(len(u_m)))]) #[u,v,w]

    b_ab = []

    for t in range(len(b[0])):
        s_dummy = np.array([s_unit_proj[0][t], s_unit_proj[1][t], s_unit_proj[2][t]])
        b_dummy = np.array([b[0][t], b[1][t], b[2][t]])
    #     print(s_dummy, b_dummy)
        b_ab.append((b_dummy.dot(s_dummy)))

    b_ab = b_ab*u.m

    V_ab = data['vis'][:,baseline,0]*u.Jy

    g_ab = V_ab/I_0*np.exp(1j*2*np.pi*b_ab*nu/c.c)
    return g_ab.decompose()

baselines = data['u_meters'].shape[1]

g_abs = []

for i in range(baselines):
    g_abs.append(g_ab(i).value)
g_abs # [baseline, time], uncomment for long list

diff_tab = Table(names=('RA', 'Diff'))

def diff(b1, b2):
    myst_data = np.load('data_mystery_fixed.npz')

    nu = myst_data['frequencies_Hz']*u.Hz
    baseline_1 = b1
    baseline_2 = b2
    t = [0,1]


    V_ab1 = myst_data['vis'][t, baseline_1, 0]
    V_ab2 = myst_data['vis'][t, baseline_2, 0]

    g_ab1 = g_abs[baseline_1][t]
    g_ab2 = g_abs[baseline_2][t]

    # myst_data['u_meters'][t, baseline_1]

    u_m1 = list(myst_data['u_meters'][t,baseline_1]) #ntimes, nbaselines
    v_m1 = list(myst_data['v_meters'][t,baseline_1]) #ntimes, nbaselines

    u_m2 = list(myst_data['u_meters'][t,baseline_2]) #ntimes, nbaselines
    v_m2 = list(myst_data['v_meters'][t,baseline_2]) #ntimes, nbaselines


    b1 = np.array([u_m1, v_m1, list(np.zeros(len(u_m1)))]) #[u,v,w]
    b2 = np.array([u_m2, v_m2, list(np.zeros(len(u_m2)))]) #[u,v,w]

    ovro = EarthLocation(lat=37.2317*u.deg,lon=-118.2951*u.deg)
    time = Time(myst_data['times_MJD'], format = 'mjd')
    aa = AltAz(obstime = time, location=ovro)

    ras = np.linspace(0, 359, 360)*u.deg

    # diff_tab = Table(names=('RA', 'Diff'))

    for ra in ras:

        qso = SkyCoord(ra = ra, dec = myst_data['DEC_degrees']*u.deg)

        qso_aa = qso.transform_to(aa)

        s = [qso_aa.az, qso_aa.alt]

        s_proj = [np.sin(s[0]), np.cos(s[0]), np.sin(s[1])]

        s_unit_proj = s_proj/np.linalg.norm(s_proj)

    #     print(s_unit_proj)

        b_ab1 = []
        b_ab2 = []

        for t in range(len(b1[0])): # could also use b2[0]
            s_dummy = np.array([s_unit_proj[0][t], s_unit_proj[1][t], s_unit_proj[2][t]])
            b1_dummy = np.array([b1[0][t], b1[1][t], b1[2][t]])
            b2_dummy = np.array([b2[0][t], b2[1][t], b2[2][t]])
        #     print(s_dummy, b_dummy)
            b_ab1.append((b1_dummy.dot(s_dummy)))
            b_ab2.append((b2_dummy.dot(s_dummy)))

        b_ab1 = b_ab1*u.m
        b_ab2 = b_ab2*u.m
    #     print(b_ab1, b_ab2)

        LHS = V_ab1/g_ab1*np.exp(1j*2*np.pi*b_ab1*nu/c.c)

        RHS = V_ab2/g_ab2*np.exp(1j*2*np.pi*b_ab2*nu/c.c)

        diff = LHS-RHS
    #     print(f'RA: {ra:.2f}, Difference: {np.abs(diff[0]):.2f}')
        diff_tab.add_row([ra, np.abs(diff[0])])

for i in range(100, 125):
    print(i)
    diff(i,i+1)
diff_tab

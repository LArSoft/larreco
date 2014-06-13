import os

cc = [ x.replace('.cc','') for x in os.listdir('.') if x.endswith(".cc")]

for s in cc:

    os.system('mv %s.cc %s.cxx' % (s,s))

hh = [x.replace('.hh','') for x in os.listdir('.') if x.endswith('.hh')]

for h in hh:

    os.system('mv %s.hh %s.h' % (h,h))

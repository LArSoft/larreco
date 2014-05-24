import os

headers=["%s" % h.replace(".cc","") for h in os.listdir(".") if h.endswith(".cc")]
for h in headers:
    os.system("mv %s.cc %s.cxx" % (h,h))

    


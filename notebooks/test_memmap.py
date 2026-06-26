import logging
from roman_cuts import RomanCuts
from roman_cuts import log as roman_log

roman_log.setLevel(logging.INFO)
roman_log.addHandler(logging.StreamHandler())

ff = ["rimtimsim_WFI_lvl02_F146_SCA02_field03_rampfitted_r1920c1920_256x256_sim.asdf"]

rcube = RomanCuts(field=FIELD, sca=SCA, filter="F146", file_list=ff, file_format="asdf")

ra, dec = 268.490404, -29.209571
rcube.make_cutout(radec=(ra, dec), size=(17, 17), dithered=False)

ra, dec = 268.49444172, -29.20351073
rcube.make_cutout(radec=(ra, dec), size=(25, 25), dithered=True)
"""Synthetic end-to-end row-layout and stage-count tests (no 6.57GB catalog)."""
import os
import tempfile
import unittest
from pathlib import Path
import numpy as np
from astropy.io import fits

from pm_footprint.local_v13_catalog_idcheck import table_specs
from pm_footprint.v13_full_catalog_paper_stages import (
    NAMES, extract_columns, selected_stages)


class TestFullCatalogStages(unittest.TestCase):
    def test_positive_nonpositive_2sigma_and_use_phot(self):
        # 1: positive F277W; 2: negative F277W (2sigma substitute)
        # 3: positive ratio EXACTLY 30 (must fail strict >30)
        # 4: eligible except use_phot false.
        n=4
        values={}
        for name in NAMES:
            if name=="ID":values[name]=np.arange(1,n+1,dtype=np.int64)
            elif name=="use_phot":values[name]=np.array([1,1,1,0],dtype=np.uint8)
            elif name.startswith("f_"):values[name]=np.zeros(n,dtype=np.float64)
            elif name.startswith(("e_","new_e_")):
                values[name]=np.ones(n,dtype=np.float64)*2
            elif name=="mag_F444W":values[name]=np.ones(n,dtype=np.float64)*25
            else:values[name]=np.ones(n,dtype=np.float64)*150
        values["f_F444W"]=np.ones(n,dtype=np.float64)*150
        values["new_e_F444W"]=np.ones(n,dtype=np.float64)*5
        values["e_F444W"]=np.ones(n,dtype=np.float64)*5
        values["f_F277W"]=np.array([2,-1,5,2],dtype=np.float64)
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/"tiny.fits"
            columns=[]
            for name in NAMES:
                data=values[name]
                fmt="K" if name=="ID" else ("B" if name=="use_phot" else "D")
                columns.append(fits.Column(name=name,format=fmt,array=data))
            fits.HDUList([fits.PrimaryHDU(),
                          fits.BinTableHDU.from_columns(columns)]).writeto(path)
            offset,rowsize,num,specs=table_specs(path)
            self.assertEqual(num,n)
            fd=os.open(path,os.O_RDONLY)
            try:
                buff=os.pread(fd,rowsize*n,offset)
            finally:
                os.close(fd)
            c=extract_columns(buff,offset,n,rowsize,specs)
            np.testing.assert_array_equal(c["ID"],np.arange(1,n+1))
            all_stages=selected_stages(c,"new_e_",False)
            filtered=selected_stages(c,"new_e_",True)
            self.assertEqual(int(all_stages["snr"].sum()),4)
            self.assertEqual(int(all_stages["ratio"].sum()),3)
            self.assertEqual(int(all_stages["mag"].sum()),3)
            self.assertEqual(int(filtered["snr"].sum()),3)
            self.assertEqual(int(filtered["ratio"].sum()),2)
            self.assertEqual(int(filtered["mag"].sum()),2)


if __name__=="__main__":unittest.main()

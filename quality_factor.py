from dataclasses import dataclass
import numpy as np
from scipy.signal import find_peaks

@dataclass
class QualityFactor:
    
    data: dict
    peak_idx: int

    @property
    def peak_wave_length(self):
        return self.data['wave length'][self.peak_idx]
    
    @property
    def peak_transmission(self):
        return self.data['transmission db'][self.peak_idx]
    
    @property
    def cut_off_point(self):
        return self.peak_transmission - 3
    
    @property
    def lh_cutoff(self):
        for idx, val in enumerate(self.data['transmission db'][0:self.peak_idx]):
            if self.data['transmission db'][self.peak_idx - idx] <= self.cut_off_point:
                return self.data['wave length'][self.peak_idx - idx]
        return None
    
    @property
    def rh_cutoff(self):
        for idx, val in enumerate(self.data['transmission db'][self.peak_idx:-1]):
            if self.data['transmission db'][self.peak_idx + idx] <= self.cut_off_point:
                return self.data['wave length'][self.peak_idx + idx]
        return None
    
    @property
    def quality_factor(self):
        if (self.lh_cutoff is not None) and (self.rh_cutoff is not None):
            return self.peak_wave_length / (self.rh_cutoff - self.lh_cutoff)
        else:
            return None

def get_quality_factor(dat: dict) -> list:
    db_peaks = find_peaks(dat['transmission db'], prominence=(10,))
    qf = []
    for idx in db_peaks[0]:
        qf.append( QualityFactor( dat, idx ) )
    return qf

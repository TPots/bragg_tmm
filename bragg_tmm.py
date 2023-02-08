from bragg import BraggGrating
import numpy as np
from numpy.linalg import matrix_power as np_mp
from dataclasses import dataclass
import quality_factor
import math
import os
import csv
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

@dataclass
class tmm:

    model:BraggGrating

    @property
    def wave_length_dependant_IOR(self) -> object:

        n_1 = self.model.FRC_1
        n_2 = self.model.FRC_2
        n_3 = self.model.FRC_3
        wl_0 = self.model.center_wave_length

        return lambda wl : n_1 + n_2 * 1e6 * (wl - wl_0) + n_3 * 1e12 * (wl - wl_0) ** 2

    @property
    def modelHash(self) -> int:
        return hash(self.model)

    def sweep(self, wl_start:float, wl_stop:float, step: float):

        step_count = round(( wl_stop - wl_start ) / step)
        span = np.linspace(wl_start, wl_stop, step_count)

        dat = {
            'wave length': [],
            'transmission amp': [],
            'reflection amp': [],
            'transmission db': [],
            'reflection db': []
        }

        for wl in span:
            t,r,t_log,r_log = self.tmm(wl)
            dat['wave length'].append(wl)
            dat['transmission amp'].append(t)
            dat['reflection amp'].append(r)
            dat['transmission db'].append(t_log)
            dat['reflection db'].append(r_log)

        return dat

    def tmm(self, wave_length: float):

        n_effective = self.wave_length_dependant_IOR(wave_length)
        n_1 = n_effective + self.model.delta_n / 2
        n_2 = n_effective - self.model.delta_n / 2

        hsm_1 = self.homogenious_section_matrix(n_1, wave_length, self.model.section_length, self.model.loss)
        hsm_2 = self.homogenious_section_matrix(n_2, wave_length, self.model.section_length, self.model.loss)
        sm_1_2 = self.step_matrix(n_1, n_2)
        sm_2_1 = self.step_matrix(n_2, n_1)
        net = np_mp(hsm_1 @ sm_1_2 @ hsm_2 @ sm_2_1, self.model.grating_count)

        hsm_fp = self.homogenious_section_matrix(n_effective, wave_length, 100e-6, 69)

        net = net@hsm_fp@net

        transmission = np.abs( 1 / net[0][0] ) ** 2
        reflection = np.abs( net[1][0] / net[0][0] ) ** 2

        return transmission, reflection, 20 * np.log10( transmission ), 20 * np.log10( reflection )

    def homogenious_section_matrix(self, n: float, wave_length: float, length:float, loss: float):

        beta = 2 * np.pi * n / wave_length - 1j * loss / 2

        mat = np.array \
        ([
            (np.exp(1j * beta * length), 0),
            (0, np.exp(-1j * beta * length))
        ])

        return mat

    def step_matrix(self, n_1: float, n_2: float):

        a = (n_1 + n_2) / (2 * np.sqrt(n_1 * n_2))
        b = (n_1 - n_2) / (2 * np.sqrt(n_1 * n_2))

        mat = np.array \
        ([
            (a, b),
            (b, a)
        ])

        return mat

    def save(self, dat:dict):
        self.save_as_csv(dat)
        self.save_as_plot(dat)
        return

    def save_as_plot(self, dat:dict):

        path = os.path.join(os.getcwd(), 'tmm')
        if not os.path.exists(path):
            os.makedirs(path)

        file_name = f'model_{self.modelHash}_{int(self.model.center_wave_length * 1e9)}_tmm.png'

        wl = dat['wave length']
        t = dat['transmission amp']
        r = dat['reflection amp']
        t_log = dat['transmission db']
        r_log = dat['reflection db']

        plt.subplot(2,1,1)
        plt.plot(wl,t,label='Transmission [amp]')
        plt.plot(wl,r,label='Reflection [amp]')
        plt.legend()

        plt.subplot(2,1,2)
        plt.plot(wl,t_log,label='Transmission [db]')
        plt.plot(wl,r_log,label='Reflection [db]')
        plt.legend()

        plt.savefig(os.path.join(path,file_name))
        print(f'Saved plot to {file_name}')
        plt.show()
        return

    def save_as_csv(self,dat:dict):
        path = os.path.join(os.getcwd(), 'tmm')
        if not os.path.exists(path):
            os.makedirs(path)

        file_name = f'model_{self.modelHash}_{int(self.model.center_wave_length * 1e9)}_tmm.csv'

        with open(os.path.join(path,file_name),'w') as dat_file:
            csv_file = csv.writer(dat_file, delimiter=',')
            csv_file.writerow(['Index','Wave Length', 'Transmission [amp]', 'Reflection [amp]', 'Transmission [db]', 'Reflection [db]'])
            for idx, wl in enumerate(dat['wave length']):
                csv_file.writerow(dat[key][idx] for key in dat)

        print(f'Saved data to {file_name}')
        return

def main():
    d_lam = 16.0346e-9
    b = BraggGrating( 1310e-9, 2.42955, -1.5789, -0.0581827 , 269e-9, 100, 133863, 441.234661241 )
    mod = tmm(b)
    dat = mod.sweep(b.center_wave_length - 2*d_lam,b.center_wave_length + 2*d_lam,1e-12)
    mod.save(dat)
    qf = quality_factor.get_quality_factor(dat)
    print(f'75%: {np.percentile([q.quality_factor for q in qf],75)}')
    return

if __name__ == '__main__':
    main()

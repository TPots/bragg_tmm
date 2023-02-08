from dataclasses import dataclass

@dataclass(frozen=True)
class BraggGrating:
    center_wave_length: float
    FRC_1 : float
    FRC_2 : float
    FRC_3 : float
    period: float
    grating_count: int
    coupling_coefficent: float
    loss: float
    
    @property
    def section_length(self) -> float:
        return self.period/2
    
    @property
    def length(self) -> float:
        return self.grating_count * self.period
    
    @property
    def delta_n(self) -> float:
        return self.coupling_coefficent * self.center_wave_length / 2
    

def main():
    b = BraggGrating(1310e-9, 2.4397, -1.84213, 0.723176 , 270e-9, 800, 98690, 3e-2)
    print(hash(b))
    return
    
    
if __name__ == '__main__':
    main()
    

CL95_results = [0.3457,0.0151,0.0093,0.0063,0.0044,0.0034,0.0024,0.0015]
numbers = [271641,296329,277381,250330,215033,187073,111094,84843]
masses = [0.005, 0.05, 0.1, 0.2, 0.5, 1, 10, 100]

print('masses = ',masses)

Nsigs=[]
for i in range(len(numbers)):
    Nsigs.append(CL95_results[i]*numbers[i])
print('Nsigs = ',Nsigs)

totalnumber = 1000000
effs=[]
for i in range(len(numbers)):
    effs.append(numbers[i] / totalnumber)
print('effs = ',effs)

NDMs=[]
for i in range(len(masses)):
    NDMs.append(3*100000 / masses[i])
print('NDMs = ',NDMs)

L=1/60.
time_year = 3*10000000
ULs=[]
for i in range(len(Nsigs)):
    ULs.append( (Nsigs[i]) / (NDMs[i]*L*time_year*effs[i]) )
print("ULs = ",ULs)

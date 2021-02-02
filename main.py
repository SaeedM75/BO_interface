from BO_interfaces.BO import *
import os
import argparse

# Command to run VASP executable.
VASP_RUN_COMMAND = 'mpirun -np 56 vasp_ncl'
# Define the name for output file.
OUTFILENAME = 'vasp.out'
# Define the path direct to the VASP pesudopotential.
VASP_PP_PATH = '/home/smoayedp/pseudopotential'


def parse_argument():
    """
    bounds: The parameter that determines the range of shift in
            Cartetian coordinates in (x,y,z) format

    kappa: The parameter to control exploration and exploitation.
           exploitation 0 <-- kappa --> 10 exploration

    threshold: Convergence threshold of Bayesian optimization process.
    """
    parser = argparse.ArgumentParser(description='params')
    parser.add_argument('--bounds', dest='bounds',
                        type=tuple, default=[2, 2, 2])
    parser.add_argument('--kappa', dest='kappa', type=float, default=4)
    parser.add_argument('--threshold', dest='threshold',
                        type=float, default=0.05)
    parser.add_argument('--import_kpath', dest='import_kpath',
                        type=bool, default=False)

    return parser.parse_args()

def main():
    args = parse_argument()
    k = args.kappa
    br = args.br
    import_kpath = args.import_kpath
    os.environ['VASP_PP_PATH'] = VASP_PP_PATH
    if os.path.exists('POSCAR_shifted'):
        os.remove('POSCAR_shifted')

    calculate(command=VASP_RUN_COMMAND, outfilename=OUTFILENAME,
              method='dftu', import_kpath=import_kpath)
    outcar = open("dftu/scf/OUTCAR")
    for l in outcar.readlines():
        if "energy without entropy" in l:
            energy = float(l.split()[7])
    with open("BO_result.txt", "w") as res_file:
        print(0, 0, 0, energy, file=res_file)

    threshold = args.threshold
    for i in range(100):
        if os.path.exists('POSCAR_shifted'):
            os.rename('POSCAR_shifted', 'POSCAR_shifted_' + str(i))

        bayesianOpt = bayesOpt_DFTU(path='./')
        obj_next = bayesianOpt.bo()
        if i == 0:
            obj_ref = obj_next
        else:
            if (np.abs(obj_next[0] - obj_ref[0]) < threshold) and (np.abs(obj_next[1] - obj_ref[1]) < threshold) and \
                    (np.abs(obj_next[2] - obj_ref[2]) < threshold):
                break
            obj_ref = obj_next

        calculate(command=VASP_RUN_COMMAND, outfilename=OUTFILENAME,
                  method='dftu', import_kpath=import_kpath)
        outcar = open("dftu/scf/OUTCAR")
        for l in outcar.readlines():
            if "energy without entropy" in l:
                energy = float(l.split()[7])
        with open("BO.txt", "a") as res_file:
            print(energy, file=res_file)


if __name__ == "__main__":
    main()

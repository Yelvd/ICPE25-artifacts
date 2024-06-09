import argparse
import math

config_string = """<?xml version="1.0" ?>
<hemocell>

<parameters>
    <warmup> 0 </warmup> <!-- Number of LBM iterations to prepare fluid field. -->
    <outputDirectory>tmp_1</outputDirectory>
    <logDirectory>log_1</logDirectory>
</parameters>

<ibm>
    <stepMaterialEvery> 20 </stepMaterialEvery> <!-- Update particle material model after this many fluid time steps. -->
    <stepParticleEvery> 5 </stepParticleEvery> <!-- Update particles position after this many fluid time steps. -->
</ibm>

<domain>
    <shearrate> 20 </shearrate>   <!--Shear rate for the fluid domain. [s^-1] [25]. -->
    <fluidEnvelope> 2 </fluidEnvelope>
    <rhoP> 1025 </rhoP>   <!--Density of the surrounding fluid, Physical units [kg/m^3]-->
    <nuP> 1.1e-6 </nuP>   <!-- Kinematic viscosity of blood plasma, physical units [m^2/s]-->
    <dx> 5.0e-7 </dx> <!--Physical length of 1 Lattice Unit -->
    <dt> -1 </dt> <!-- Time step for the LBM system. A negative value will set Tau=1 and calc. the corresponding time-step. -->
    <refDir> 1 </refDir>   <!-- Used for resloution  setting and  Re calculation as well -->
    <nx> {x} </nx>  <!-- Number of numerical cell in the reference direction -->
    <ny> {y} </ny>  <!-- Number of numerical cell in the reference direction -->
    <nz> {z} </nz>  <!-- Number of numerical cell in the reference direction -->
    <blockSize> -1 </blockSize>
    <kBT> 4.100531391e-21 </kBT> <!-- in SI, m2 kg s-2 (or J) for T=300 -->
    <particleEnvelope> 25 </particleEnvelope>
</domain>

<sim>
    <tmax> {iters} </tmax> <!-- total number of iterations -->
    <tmeas> 100000 </tmeas> <!-- interval after which data is written -->
</sim>

<benchmark>
    <tslice> 500000 </tslice>
    <imbalance> 0 </imbalance>
    <writeOutput> 0 </writeOutput>
</benchmark>

</hemocell>
"""

def create_primes(n):
    primes = []
    # Print the number of two's that divide n
    while n % 2 == 0:
        primes.append(2)
        n = n // 2
         
    # n must be odd at this point
    # so a skip of 2 ( i = i + 2) can be used
    for i in range(3,int(math.sqrt(n))+1,2):
         
        # while i divides n , print i ad divide n
        while n % i== 0:
            primes.append(i)
            n = n // i
             
    # Condition if n is a prime
    # number greater than 2
    if n > 2:
        primes.append(n)
        print(n)

    return primes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xs", type=int, help="x size of atomic-block")
    parser.add_argument("-y", "--ys", type=int, help="y size of atomic-block")
    parser.add_argument("-z", "--zs", type=int, help="z size of atomic-block")
    parser.add_argument("-n", "--tasks", type=int, help="number of tasks")
    parser.add_argument("-i", "--iters", type=int, help="number of tasks", default=10000)
    parser.add_argument("-f", "--out_file", type=str, help="name of the output file", default=None)
    args = parser.parse_args()

    if args.out_file is None:
        args.out_file = f"{args.xs}-{args.ys}-{args.zs}_{args.tasks}_config.xml"

    primes = create_primes(args.xs*args.ys*args.zs*args.tasks)
    primes.sort(reverse=True)

    dims = [1, 1, 1]

    for i, p in enumerate(primes):
        dims[i % 3] = dims[i % 3] * p
        
    dims.sort(reverse=True)
    print(dims)


    print(args.out_file)

    f = open(args.out_file, "w")
    f.write(config_string.format(x=dims[0], y=dims[1], z=dims[2], iters=args.iters))
    f.close()


if __name__ == "__main__":
    main()

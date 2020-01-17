#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from mykit.core.symmetry import primitize, is_primitive
from mykit.vasp.poscar import Poscar


def pv_primitize():
    """Primizite the VASP POSCAR structure
    """
    parser = ArgumentParser(
        description=pv_primitize.__doc__, formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument("-f", default="POSCAR", help="POSCAR input file")
    parser.add_argument("-o", dest="output", default="POSCAR_prim", help="output path")
    parser.add_argument("-D", dest="debug", action="store_true", help="debug mode")
    args = parser.parse_args()

    c = Poscar.read_from_file(args.f)
    if args.debug:
        print(c)
        print(type(c))
    flag = is_primitive(c)

    if flag:
        print("%s is already a primitive cell." % args.f)
        return
    elif flag is None:
        print("%s is a bad input." % args.f)
        return
    
    c = primitize(c)
    print("Primitizing %s and write to %s" % (args.f, args.output))
    c.write(args.output)


if __name__ == "__main__":
    pv_primitize()

import cli
import sys

def main():
    # Parse arguments (from command line or programmatic input)
    args = cli.parse_args()

    if args.method == "blast":
        try:
            import blast_101_search as blast
            blast.blast101_run()
            sys.exit(0)
        except Exception as e:
            # Log the error message and exit with a non-zero exit code.
            print(f"Error running BLAST pipeline: {e}", file=sys.stderr )
            sys.exit(1)
    elif args.method == "sw":
        try:
            import smith_waterman_search
            smith_waterman_search.smith_run()
            sys.exit(0)
        except Exception as e:
            print(f"Error running Smith-Waterman pipeline: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.method == "stats":
        try:
            import calc_bit_and_evalues
            calc_bit_and_evalues.build_fit()
            sys.exit(0)
        except Exception as e:
            print(f"Error running statistical analysis: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.method == "tests":

        import unittest
        print("**********************************************************************************")
        print("*                                                                                *")
        print("*                            Running Tests                                       *")
        print("*                                                                                *")
        print("**********************************************************************************\n")

        suite = unittest.TestLoader().discover('.', pattern='test_*.py')

        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)

        sys.exit(not result.wasSuccessful())

if __name__ == "__main__":
    main()
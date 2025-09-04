#     Alexander invariants for knots from PD codes.
#     This module computes generalizations of the classical Alexander polynomial,
#     as introduced by Alexander. It includes the Alexander polynomials, the Alexander
#     Invariants, and the primary decomposition.

#     Usage:
#       In SageMath, run:
#         load('alex.sage')
#     Then call:
#        alexander_invariants(pd_code)
#
#     Example:
#         pd_code = [[1,5,2,4],[3,1,4,6],[5,3,6,2]]
#         alexander_invariants(pd_code, decomp='alex_polys_decomp')  (Default parameter value.)
#         alexander_invariants(pd_code, decomp='primary_decomp')
#         alexander_invariants(pd_code, decomp='invariant_factor_decomp')


import itertools
R_ZZ = PolynomialRing(ZZ, 't')
t = R_ZZ.gen()
# --- decomposition mode constants ---
PRIMARY_DECOMP = "primary_decomp"
INVARIANT_FACTOR_DECOMP = "invariant_factor_decomp"
ALEX_POLYS_DECOMP = "alex_polys_decomp"

def _check_pd_code(pd_code): 
#    Check if a PD code is valid.
#    A PD code is valid if:
#    - Each entry is a list of 4 integers [a, b, c, d]
#    - c - a ≡ 1 mod 2N
#    - d - b ≡ ±1 mod 2N
#    - Each label appears exactly twice
#
#    INPUT:
#    - ``pd_code`` -- list of 4-tuples
#
#    OUTPUT:
#    - ``True`` if valid; raises ``ValueError if invalid`` otherwise

    for x in pd_code:
        if len(x) != 4:
            raise ValueError("Invalid PD code")
    N = 2*len(pd_code)
    if not all(
        ( (x[2] - x[0]) % N == 1 and (x[3] - x[1]) % N in (1, N - 1) )
        for x in pd_code):
        raise ValueError("Invalid PD code")
    
    counts = [0]*N
    for pdc in pd_code:
        for x in pdc:
            counts[x % N] += 1
    if not all(count == 2 for count in counts):
        raise ValueError("Invalid PD code")

    return True
    

def alexander_invariants(pd_code, decomp = 'alex_polys_decomp'):
#    Returns Alexander invariants of a knot, in any of the forms 'primary 
#    decomposition', 'alexander invariants decomposition', or 'alexander 
#    polynomial decomposition'.
#
#    EXAMPLES::
#        sage: pd_code = [[6,2,7,1],[8,3,9,4],[16,11,1,12],[2,14,3,13],
#        ....:            [4,15,5,16],[10,6,11,5],[12,7,13,8],[14,10,15,9]]
#        sage: alexander_invariants(pd_code, decomp='primary_decomp')
#        [[t^2 - 3*t + 1, [1]], [t^2 - t + 1, [1, 1]]]
#
#        sage: alexander_invariants(pd_code, decomp='invariant_factor_decomp')
#        [(t^2 - 3*t + 1) * (t^2 - t + 1), t^2 - t + 1]
#
#        sage: alexander_invariants(pd_code, decomp='alex_polys_decomp')
#        [(t^2 - 3*t + 1) * (t^2 - t + 1)^2, t^2 - t + 1]
#
#        sage: bad_pd_code = [[6,2,8,1],[8,3,9,4],[16,11,1,12],[2,14,3,13],
#        ....:                [4,15,5,16],[10,6,11,5],[12,7,13,8],[14,10,15,9]]
#        sage: try:
#        ....:     alexander_invariants(bad_pd_code, decomp='primary_decomp') #raises ValueError
#        ....: except ValueError as e:
#        ....:     print(e)
#        Invalid PD code
#
#        sage: _check_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]])
#        True
#
#        sage: _check_pd_code([[0, 2, 1]])
#        Traceback (most recent call last):
#        ...
#        ValueError: Invalid PD code

    _check_pd_code(pd_code) #will raise error if invalid
    return AlexanderInvariants.alex_polys(pd_code, decomp)

def flatten_AI(nested_list):
    return list(itertools.chain.from_iterable(nested_list))


class AlexanderInvariants:
    @staticmethod
    def alex_polys(pd, decomp='alex_polys_decomp'):
        """
        Compute Alexander polynomials of a knot from its PD code.

        Parameters:
            pd : list
                Planar diagram data describing the knot.
            decomp : str, optional
                Type of decomposition to return:
                'primary_decomp' (default), 'invariant_factor_decomp', or 'alex_polys_decomp'.

        Returns:
            List of polynomials according to decomposition type.
            'primary_decomp' = primary decomposition of Alexander module
            'invariant_factor_decomp' = invariant factors of Alexander module
            'alex_polys_decomp' = list of Alexander polynomials with exponents of summands
        """
        fm = AlexanderInvariants._pd_to_fox_matrix(pd)   #computes fox_matrix from PD_code
        alex_poly = AlexanderInvariants._alex_poly_from_fm(fm)  #computes Alexander polynomial
        alex_polys = AlexanderInvariants._find_possible_primary_decomps_all_factors(fm)  
        #find splitting, expressed in three ways
        if decomp == PRIMARY_DECOMP:
            return alex_polys[0]
        if decomp == INVARIANT_FACTOR_DECOMP:
            return [x.factor() for x in alex_polys[1]]
        if decomp == ALEX_POLYS_DECOMP:
            return alex_polys[2]

##########CODE TO COMPUTE FOX MATRIX FROM PD_CODE############
    @staticmethod
    def _pd_to_fox_matrix(pd):
        """
        Convert PD code to the Fox matrix used for Alexander polynomial calculation.

        Parameters:
            pd : list
                Planar diagram data.

        Returns:
            Matrix over ZZ[t] representing the Fox calculus relations.
        """
        cross_num = len(pd)

        def _normalize_pd(pd):
            """Normalize all arc indices modulo 2 * number of crossings."""
            return [[y % (2  * cross_num) for y in x] for x in pd]

        def _is_right_handed(cr_data):
            """
            Determine if a crossing is right-handed based on arc ordering.

            Returns:
                bool: True if right-handed, False otherwise.
            """
            return (cr_data[3] - cr_data[1] == 1 ) or (cr_data[3] - cr_data[1] < -1 )

        def _pre_wirtinger_from_cross(cr_data):
            """
            Produce a preliminary Wirtinger relation from crossing data,
            adjusting based on crossing handedness.
            """
            if _is_right_handed(cr_data):
                return [cr_data[1], cr_data[2], cr_data[3], cr_data[0]]
            else:
                return [cr_data[1], cr_data[0], cr_data[3], cr_data[2]]

        def _build_pre_wirtinger_relations(pd_normalized):
            """Build all preliminary Wirtinger relations for normalized PD code."""
            return [_pre_wirtinger_from_cross(cr_data) for cr_data in pd_normalized]

        def _unify_over_arcs(pre_wirt):
            """
            Unify arcs appearing as 'over' arcs to create consistent labeling,
            necessary for Fox matrix construction.
            """
            for i in range(cross_num):
                in_over = pre_wirt[i][0]
                out_over = pre_wirt[i][2]
                # Replace all occurrences of out_over by in_over to unify arcs
                pre_wirt = [[in_over if x == out_over else x for x in row] for row in pre_wirt]
            return pre_wirt

        def _minimize_wirt(wirt_relations):
            """
            Renumber arcs to a minimal consecutive index set for matrix indexing.
            """
            values = list(set(flatten_AI(wirt_relations)))
            index_map = {val: idx for idx, val in enumerate(values)}
            return [[index_map[x] for x in rel] for rel in wirt_relations]

        def _build_fox_matrix(min_wirt):
            """
            Construct the Fox matrix from minimized Wirtinger relations.
            Matrix entries are polynomials in ZZ[t].
            """
            fm = zero_matrix(R_ZZ, cross_num - 1, cross_num - 1 )
            for i in range(cross_num - 1 ):
                rel = min_wirt[i]
                # Assign entries according to Fox calculus rules
                if rel[0] != 0 :
                    fm[i, rel[0 ] - 1] = 1  - t
                if rel[1] != 0 :
                    fm[i, rel[1 ] - 1 ] = t
                if rel[3] != 0 :
                    fm[i, rel[3] - 1 ] = -1 
            return fm

        # --- Execution flow for conversion ---
        pd_normalized = _normalize_pd(pd)
        pre_wirt = _build_pre_wirtinger_relations(pd_normalized)
        pre_wirt = _unify_over_arcs(pre_wirt)
        min_wirt = _minimize_wirt(pre_wirt)
        fox_mat = _build_fox_matrix(min_wirt)

        return fox_mat

##########CODE TO COMPUTE ALEX POLY AND ITS FACTORIZATION NORMALIZED############
    @staticmethod
    def _alex_poly_from_fm(fm):
        """
        Compute the Alexander polynomial by taking the determinant
        of the Fox matrix and normalising it.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].

        Returns:
            Polynomial in ZZ[t] representing the Alexander polynomial.
        """
        alex_poly = AlexanderInvariants._reduce_poly(det(fm))  #"reduce_poly" normalizes polynomial
        return R_ZZ(alex_poly)


    @staticmethod
    def _alex_poly_factors_from_fm(fm):
        """
        Factor the Alexander polynomial obtained from Fox matrix.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].

        Returns:
            List of (factor, exponent) pairs.
        """
        alex = AlexanderInvariants._alex_poly_from_fm(fm)
        factors = alex.factor()
        return list(factors)

    @staticmethod
    def _reduce_poly(f):
        """
        Normalize polynomial f by dividing out all factors of t,
        then make the constant term positive.

        Parameters:
            f : Polynomial
                Polynomial in ZZ[t].

        Returns:
            Reduced polynomial with positive constant term.
        """
        t = R_ZZ.gen()
        g = R_ZZ(f)

        t_order = g.valuation(t)
        if t_order > 0 :
            g = g // (t**t_order)

        return g if g.subs(t=0 ) > 0  else -g


##########CODE TO FIND ALL POSSIBLE PRIMARY DECOMPOSITIONS############
    @staticmethod
    def _find_possible_primary_decomp_factor(fm, factor):
        """
        Analyze possible splittings for a given factor of the Alexander polynomial.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].
            factor : tuple
                Polynomial factor and its exponent.

        Returns:
            List of possible exponent splittings.
        """
        poly = factor[0 ]
        exp = factor[1]

        if exp == 1:
            out_value = [poly, [[1]]]
        else:
            poss_split_based_on_num_field = AlexanderInvariants._exponent_gt_one_num_theory(fm, factor)
            if len(poss_split_based_on_num_field[1]) == 1:
                out_value = poss_split_based_on_num_field
            else:
                coker_dim = AlexanderInvariants._cokernel_degree(fm, poly)
                current_poss = poss_split_based_on_num_field[1]
                # Filter exponent splittings to those matching cokernel dimension
                new_poss = [x for x in current_poss if max(x) == coker_dim]
                out_value = [factor[0], new_poss]

        return out_value


    @staticmethod
    def _find_possible_primary_decomps_all_factors(fm):
        """
        Find all possible splittings of the Alexander polynomial
        based on factorization and algebraic decomposition.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].

        Returns:
            List containing splitting info, primary decomposition, and Alexander polys.
        """
        alex_poly = AlexanderInvariants._alex_poly_from_fm(fm)

        if alex_poly == 1:
            # Properly return consistent formats for all three types
            primary_decomp = [[R_ZZ(1), [1]]]
            alex_invariants = [R_ZZ(1)]
            alex_polys = [R_ZZ(1)]
            return [primary_decomp, alex_invariants, alex_polys]

        primary_decomp = AlexanderInvariants._get_primary_splitting(fm)
        if primary_decomp == "underdetermined":
            return ["underdetermined", None, None]

        alex_invariants = AlexanderInvariants._format_primary_decomposition(primary_decomp)
        alex_polys = AlexanderInvariants._compute_alexander_polys_from_primary(alex_invariants)

        return [primary_decomp, alex_invariants, alex_polys]


    @staticmethod
    def _get_primary_splitting(fm):
        """
        Analyze the factors of the Alexander polynomial and
        determine the possible splittings.
        This is the primary decomposition.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].

        Returns:
            List of factors and their exponents or "underdetermined".
            Each element in list is of form [f(t), [[3,3,1],[3,2,2]]] 
        """
        factors = AlexanderInvariants._alex_poly_factors_from_fm(fm)
        #factors will be a list of irred. factors and exponents.
        if not factors:
            return [[R_ZZ(1), [1]]]

        all_splittings = [
            AlexanderInvariants._find_possible_primary_decomp_factor(fm, factor)
            for factor in factors
        ]

        num_exps_possible = [len(x[1]) for x in all_splittings]
        
        #if splitting is determined, returns list with entries of form [f(t), [2,2,1]]
        #
        if AlexanderInvariants._max_or_zero(num_exps_possible) == 1:
            return [[f[0], f[1][0]] for f in all_splittings]
        else:
            return "underdetermined"

    @staticmethod
    def _compute_alexander_polys_from_primary(primary_decomp):
        """
        Compute Alexander polynomials from the primary decomposition.

        Parameters:
            primary_decomp : list
                List of polynomials representing the primary decomposition.

        Returns:
            List of Alexander polynomials.
        """
        return [prod(primary_decomp[i:]).factor() for i in range(len(primary_decomp))]

    @staticmethod
    def _format_primary_decomposition(splitting_info):
        """
        Format the primary decomposition into a list of polynomial factors raised to powers.

        Parameters:
            splitting_info : list
                List of tuples with factors and exponent lists.

        Returns:
            List of polynomials corresponding to the decomposition.
        """
        if splitting_info == [1]:
            return [1, [1]]

        max_len = max(len(x[1]) for x in splitting_info)

        padded = []
        for factor, exps in splitting_info:
            # Pad exponent lists to uniform length with zeros
            padded_exps = exps + [0] * (max_len - len(exps))
            padded.append([factor, padded_exps])

        pfactors = [x[0] for x in padded]
        zipped_exponents = list(zip(*[x[1] for x in padded]))
        #if padded = [[f, [2,1]],[g,[5,2]] then zipped_exponents = [[2,5],[1,2]]
        
        # Construct polynomials by raising factors to exponents and multiplying
        return [prod(b ** e for b, e in zip(pfactors, exp)) for exp in zipped_exponents]

    @staticmethod
    def _max_or_zero(li):
        """
        Return max of list or 0 if empty.

        Parameters:
            li : list
                List of comparable values.

        Returns:
            Maximum value or zero.
        """
        if li:
            return max(li)
        else:
            return 0 

    @staticmethod
    def _exponent_gt_one_num_theory(fm, factor):
        """
        Use number theory on a number field extension defined by the factor
        to analyze possible exponent splittings for factors with exponent > 1.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].
            factor : tuple
                Factor polynomial and exponent.

        Returns:
            List with polynomial and list of possible exponent compositions.
        """
        factor_poly = factor[0]
        factor_exp = factor[1]

        factor_poly_Z = R_ZZ(factor_poly)
        fm_Z = Matrix(R_ZZ, fm)

        # Define number field by adjoining root of factor_poly
        K = NumberField(factor_poly_Z, names=('a',)); (a,) = K._first_ngens(1)

        # Ring homomorphism from ZZ[t] to number field
        phi = R_ZZ.hom([a], K)
        fm_K = fm_Z.apply_map(phi)

        rank_fm_K = fm_K.rank()
        dim_fm = fm.nrows()
        dim_factor_part = dim_fm - rank_fm_K

        # Generate all compositions of the exponent with length equal to factor dimension
        compositions = Compositions(factor_exp, length=dim_factor_part)

        # Filter to non-increasing compositions (descending order)
        non_decreasing = [
            list(comp)
            for comp in compositions
            if sorted(comp, reverse=True) == list(comp)
        ]

        return [factor_poly, non_decreasing]

    @staticmethod
    def _cokernel_degree(fm, factor):
        """
        Compute the degree of the cokernel corresponding to a factor.

        Parameters:
            fm : Matrix
                Fox matrix over ZZ[t].
            factor : Polynomial
                Polynomial factor of Alexander polynomial.

        Returns:
            Integer valuation of the cokernel degree.
        """
        try:
            fm_inv = fm.inverse()
        except ZeroDivisionError:
            raise ValueError("Fox matrix is not invertible — likely invalid PD code.")

        lcm_denominators = lcm([
            entry.denominator()
            for row in fm_inv.rows()
            for entry in row
        ])

        e = lcm_denominators.valuation(factor)
        return e


############HELPERS TO BUILD AND CALL SQLITE DATABASE3#######
#import sqlite3
#import ast
#import os
#
#    # ---- File paths ----
#    names_file = "names.txt"
#    pd_file = "pd.txt"
#    db_file = "pd_database.sqlite"
#    
#    # ---- Create database ----
#    if os.path.exists(db_file):
#        os.remove(db_file)  # Remove old database for a clean rebuild
#    
#    conn = sqlite3.connect(db_file)
#    cursor = conn.cursor()
#    
#    # ---- Create table ----
#    cursor.execute("""
#        CREATE TABLE pd_data (
#            id INTEGER PRIMARY KEY,
#            name TEXT UNIQUE,
#            pd TEXT
#        )
#    """)
#    
#    cursor.execute("CREATE INDEX idx_name ON pd_data(name)")
#    
#    # ---- Read and insert data ----
#    with open(names_file, 'r') as nf, open(pd_file, 'r') as pf:
#        for idx, (name_line, pd_line) in enumerate(zip(nf, pf), start=1):  # <-- Start from 1
#            name = name_line.strip()
#            pd = pd_line.strip()
#            cursor.execute("INSERT INTO pd_data (id, name, pd) VALUES (?, ?, ?)", (idx, name, pd))
#    
#    conn.commit()

# ---- Define search functions ----
import sqlite3
import ast
import os

db_file = "pd_database.sqlite"  # your database file path

# Immediate check when the file is loaded
if not os.path.exists(db_file):
    print(f"""
[Warning] The database file "{db_file}" was not found.

Some functionality (the "find" functions) of this program will not work without it.

To fix this:
  1. Download and unpack the database from the github site:
     "github.com/chuck-livington/alexander-invariants".
  2. Place it in the same directory as this program, 
     or update 'db_file' to point to its location.
""")

def get_connection():
    """Safely return a new SQLite connection if database file exists."""
    if not os.path.exists(db_file):
        raise FileNotFoundError(f"Database file {db_file} not found")
    return sqlite3.connect(db_file)

#conn = sqlite3.connect(db_file)

def find_name_from_index(n):
    n = int(n)  # Convert Sage Integer to Python int
    with get_connection() as conn:  # ensures connection closes automatically
        cur = conn.execute("SELECT name FROM pd_data WHERE id = ?", (n,))
        result = cur.fetchone()
    return result[0] if result else None

def find_pd_from_index(n):
    n = int(n)
    with get_connection() as conn:
        cur = conn.execute("SELECT pd FROM pd_data WHERE id = ?", (n,))
        result = cur.fetchone()
    if result:
        return ast.literal_eval(result[0])  # convert string to Python list
    return None

def find_index_from_name(name):
    with get_connection() as conn:
        cur = conn.execute("SELECT id FROM pd_data WHERE name = ?", (name,))
        result = cur.fetchone()
    return result[0] if result else None

def find_pd_from_name(name):
    ind = find_index_from_name(name)
    if ind is not None:
        return find_pd_from_index(ind)
    return None
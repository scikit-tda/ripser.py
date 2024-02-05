/*

Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

MIT License

Original Copyright 2015-2018 Ulrich Bauer.
Modifications Copyright 2018 Christopher Tralie


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the author of this software, without
imposing a separate written license agreement for such Enhancements, then you
hereby grant the following license: a non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.


*/

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>

// #define INDICATE_PROGRESS

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
static const std::string clear_line("\r\033[K");
#endif

#if defined(USE_ROBINHOOD_HASHMAP)
#include <robin_hood.h>

template <class Key, class T, class H, class E>
using hash_map_imp = robin_hood::unordered_map<Key, T, H, E>;
#else
template <class Key, class T, class H, class E>
using hash_map_imp = std::unordered_map<Key, T, H, E>;
#endif

template <class Key, class T, class H, class E>
class hash_map : public hash_map_imp<Key, T, H, E>
{
};

typedef float value_t;
typedef int64_t index_t;
typedef int16_t coefficient_t;

static const size_t num_coefficient_bits = 8;

// 1L on windows is ALWAYS 32 bits, when on unix systems is pointer size
static const index_t max_simplex_index =
    (uintptr_t(1) << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

void check_overflow(index_t i)
{
    if
#ifdef USE_COEFFICIENTS
        (i > max_simplex_index)
#else
        (i < 0)
#endif
        throw std::overflow_error(
            "simplex index " + std::to_string((uint64_t) i) +
            " in filtration is larger than maximum index " +
            std::to_string(max_simplex_index));
}

class binomial_coeff_table
{
    /* Using flatten matrix */
    std::vector<index_t> B;
    size_t offset;

public:
    binomial_coeff_table(index_t n, index_t k) : B((n + 1) * (k + 1))
    {
        offset = k + 1;
        for (index_t i = 0; i <= n; ++i) {
            B[i * offset] = 1;
            for (index_t j = 1; j < std::min(i, k + 1); ++j)
                B[i * offset + j] =
                    B[(i - 1) * offset + j - 1] + B[(i - 1) * offset + j];
            if (i <= k)
                B[i * offset + i] = 1;
            check_overflow(B[i * offset + std::min(i >> 1, k)]);
        }
    }

    index_t operator()(index_t n, index_t k) const
    {
        assert(n < (B.size() / offset) && k < offset && n >= k - 1);
        return B[n * offset + k];
    }
};

/* Modulo operator is expensive, using a mask when modulus is equal 2
 * is much less expesive and speed-ups where observed
 */
const coefficient_t get_modulo(const coefficient_t val,
                               const coefficient_t modulus)
{
    return (modulus == 2) ? val & 1 : val % modulus;
}

coefficient_t normalize(const coefficient_t n, const coefficient_t modulus)
{
    return n > modulus / 2 ? n - modulus : n;
}

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m)
{
    std::vector<coefficient_t> inverse(m);
    inverse[1] = 1;
    // m = a * (m / a) + m % a
    // Multipying with inverse(a) * inverse(m % a):
    // 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
    for (coefficient_t a = 2; a < m; ++a)
        inverse[a] = m - (inverse[m % a] * (m / a)) % m;
    return inverse;
}

#ifdef USE_COEFFICIENTS

// https://stackoverflow.com/a/3312896/13339777
#ifdef _MSC_VER
#define PACK( ... ) __pragma( pack(push, 1) ) __VA_ARGS__ __pragma( pack(pop))
#else
#define PACK( ... ) __attribute__((__packed__)) __VA_ARGS__
#endif

PACK(struct entry_t {
    index_t index : 8 * sizeof(index_t) - num_coefficient_bits;
    index_t coefficient : num_coefficient_bits;
    entry_t(index_t _index, coefficient_t _coefficient)
        : index(_index), coefficient(_coefficient)
    {
    }
    entry_t(index_t _index) : index(_index), coefficient(0) {}
    entry_t() : index(0), coefficient(0) {}
});

static_assert(sizeof(entry_t) == sizeof(index_t),
              "size of entry_t is not the same as index_t");

entry_t make_entry(index_t i, coefficient_t c) { return entry_t(i, c); }
index_t get_index(const entry_t& e) { return e.index; }
index_t get_coefficient(const entry_t& e) { return e.coefficient; }
void set_coefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

std::ostream& operator<<(std::ostream& stream, const entry_t& e)
{
    stream << get_index(e) << ":" << get_coefficient(e);
    return stream;
}

#else

typedef index_t entry_t;
const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value)
{
    return entry_t(_index);
}
void set_coefficient(entry_t& e, const coefficient_t c) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }

typedef std::pair<value_t, index_t> diameter_index_t;
value_t get_diameter(const diameter_index_t& i) { return i.first; }
index_t get_index(const diameter_index_t& i) { return i.second; }

typedef std::pair<index_t, value_t> index_diameter_t;
index_t get_index(const index_diameter_t& i) { return i.first; }
value_t get_diameter(const index_diameter_t& i) { return i.second; }

struct diameter_entry_t : std::pair<value_t, entry_t> {
    using std::pair<value_t, entry_t>::pair;
    diameter_entry_t() {}
    diameter_entry_t(value_t _diameter, index_t _index,
                     coefficient_t _coefficient)
        : diameter_entry_t(_diameter, make_entry(_index, _coefficient))
    {
    }
    diameter_entry_t(const diameter_index_t& _diameter_index,
                     coefficient_t _coefficient)
        : diameter_entry_t(get_diameter(_diameter_index),
                           make_entry(get_index(_diameter_index), _coefficient))
    {
    }
    diameter_entry_t(const index_t& _index) : diameter_entry_t(0, _index, 0) {}
};

const entry_t& get_entry(const diameter_entry_t& p) { return p.second; }
entry_t& get_entry(diameter_entry_t& p) { return p.second; }
const index_t get_index(const diameter_entry_t& p)
{
    return get_index(get_entry(p));
}
const coefficient_t get_coefficient(const diameter_entry_t& p)
{
    return get_coefficient(get_entry(p));
}
const value_t& get_diameter(const diameter_entry_t& p) { return p.first; }
void set_coefficient(diameter_entry_t& p, const coefficient_t c)
{
    set_coefficient(get_entry(p), c);
}

template <typename Entry>
struct greater_diameter_or_smaller_index {
    bool operator()(const Entry& a, const Entry& b)
    {
        return (get_diameter(a) > get_diameter(b)) ||
               ((get_diameter(a) == get_diameter(b)) &&
                (get_index(a) < get_index(b)));
    }
};

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout>
class compressed_distance_matrix
{
public:
    std::vector<value_t> distances;
    std::vector<value_t*> rows;

    compressed_distance_matrix(std::vector<value_t>&& _distances)
        : distances(std::move(_distances)),
          rows((1 + std::sqrt(1 + 8 * distances.size())) / 2)
    {
        assert(distances.size() == size() * (size() - 1) / 2);
        init_rows();
    }

    template <typename DistanceMatrix>
    compressed_distance_matrix(const DistanceMatrix& mat)
        : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size())
    {
        init_rows();

        for (size_t i = 1; i < size(); ++i)
            for (size_t j = 0; j < i; ++j)
                rows[i][j] = mat(i, j);
    }

    value_t operator()(const index_t i, const index_t j) const;

    size_t size() const { return rows.size(); }
    void init_rows();
};

typedef compressed_distance_matrix<LOWER_TRIANGULAR>
    compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR>
    compressed_upper_distance_matrix;

template <>
void compressed_lower_distance_matrix::init_rows()
{
    value_t* pointer = &distances[0];
    for (size_t i = 1; i < size(); ++i) {
        rows[i] = pointer;
        pointer += i;
    }
}

template <>
void compressed_upper_distance_matrix::init_rows()
{
    value_t* pointer = &distances[0] - 1;
    for (size_t i = 0; i < size() - 1; ++i) {
        rows[i] = pointer;
        pointer += size() - i - 2;
    }
}

template <>
value_t compressed_lower_distance_matrix::operator()(const index_t i,
                                                     const index_t j) const
{
    return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

template <>
value_t compressed_upper_distance_matrix::operator()(const index_t i,
                                                     const index_t j) const
{
    return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}

struct sparse_distance_matrix {
    std::vector<std::vector<index_diameter_t>> neighbors;
    std::vector<value_t> vertex_births;
    index_t num_edges;

    mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator>
        neighbor_it;
    mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator>
        neighbor_end;

    sparse_distance_matrix(
        std::vector<std::vector<index_diameter_t>>&& _neighbors,
        index_t _num_edges)
        : neighbors(std::move(_neighbors)), vertex_births(_neighbors.size(), 0),
          num_edges(_num_edges)
    {
    }

    template <typename DistanceMatrix>
    sparse_distance_matrix(const DistanceMatrix& mat, const value_t threshold)
        : neighbors(mat.size()), vertex_births(mat.size(), 0), num_edges(0)
    {
        for (size_t i = 0; i < size(); ++i)
            for (size_t j = 0; j < size(); ++j)
                if (i != j && mat(i, j) <= threshold) {
                    ++num_edges;
                    neighbors[i].push_back({j, mat(i, j)});
                }
    }

    // Initialize from COO format
    sparse_distance_matrix(int* I, int* J, value_t* V, int NEdges, int N,
                           const value_t threshold)
        : neighbors(N), vertex_births(N, 0), num_edges(0)
    {
        int i, j;
        value_t val;
        for (int idx = 0; idx < NEdges; idx++) {
            i = I[idx];
            j = J[idx];
            val = V[idx];
            if (i < j && val <= threshold) {
                neighbors[i].push_back(std::make_pair(j, val));
                neighbors[j].push_back(std::make_pair(i, val));
                ++num_edges;
            } else if (i == j) {
                vertex_births[i] = val;
            }
        }

        for (size_t i = 0; i < neighbors.size(); ++i)
            std::sort(neighbors[i].begin(), neighbors[i].end());
    }

    size_t size() const { return neighbors.size(); }
};

class union_find
{
    std::vector<index_t> parent;
    std::vector<uint8_t> rank;
    std::vector<value_t> birth;

public:
    union_find(index_t n) : parent(n), rank(n, 0), birth(n, 0)
    {
        for (index_t i = 0; i < n; ++i)
            parent[i] = i;
    }

    void set_birth(index_t i, value_t val) { birth[i] = val; }

    value_t get_birth(index_t i) { return birth[i]; }

    index_t find(index_t x)
    {
        index_t y = x, z = parent[y];
        while (z != y) {
            y = z;
            z = parent[y];
        }
        y = parent[x];
        while (z != y) {
            parent[x] = z;
            x = y;
            y = parent[x];
        }
        return z;
    }

    void link(index_t x, index_t y)
    {
        x = find(x);
        y = find(y);
        if (x == y)
            return;
        if (rank[x] > rank[y]) {
            parent[y] = x;
            birth[x] = std::min(birth[x], birth[y]);  // Elder rule
        } else {
            parent[x] = y;
            birth[y] = std::min(birth[x], birth[y]);  // Elder rule
            if (rank[x] == rank[y])
                ++rank[y];
        }
    }
};

template <typename T>
T begin(std::pair<T, T>& p)
{
    return p.first;
}
template <typename T>
T end(std::pair<T, T>& p)
{
    return p.second;
}

template <typename ValueType>
class compressed_sparse_matrix
{
    std::vector<size_t> bounds;
    std::vector<ValueType> entries;

    typedef typename std::vector<ValueType>::iterator iterator;
    typedef std::pair<iterator, iterator> iterator_pair;

public:
    size_t size() const { return bounds.size(); }

    iterator_pair subrange(const index_t index)
    {
        return {entries.begin() + (index == 0 ? 0 : bounds[index - 1]),
                entries.begin() + bounds[index]};
    }

    void append_column() { bounds.push_back(entries.size()); }

    void push_back(const ValueType e)
    {
        assert(0 < size());
        entries.push_back(e);
        ++bounds.back();
    }

    void pop_back()
    {
        assert(0 < size());
        entries.pop_back();
        --bounds.back();
    }
};

/* This is the data structure from which the results of running ripser can be
 * returned */
typedef struct {
    /* The first variable is a vector of unrolled persistence diagrams
       so, for example births_and_deaths_by_dim[0] contains a list of
                [birth0, death0, birth1, death1, ..., birthk, deathk]
       for k points in the 0D persistence diagram
       and likewise for d-dimensional persistence in births_and_deaths_by_dim[d]
    */
    std::vector<std::vector<value_t>> births_and_deaths_by_dim;
    /*
      The second variable is a vector of representative cocycles for each
      dimension. For now, only cocycles above dimension 0 are added, so
      dimension 0 is an empty list For the others, cocycles_by_dim[d] holds an
      array of representative cocycles for dimension d which is parallel with
      the array of births/deaths for dimension d. Each element of the array is
      itself an array of unrolled information about the cocycle For dimension 1,
      for example, the zeroeth element of the array contains [ccl0_simplex0_idx0
      ccl0_simplex0_idx1 ccl0_simplex0_val, ccl0_simplex1_idx0
      ccl0_simplex1_idx1 ccl0_simplex1_val, ... ccl0_simplexk_idx0
      ccl0_simplexk_idx1 ccl0_simplexk_val] for a cocycle representing the first
      persistence point, which has k simplices with nonzero values in the
      representative cocycle
    */
    std::vector<std::vector<std::vector<int>>> cocycles_by_dim;
    /* The third variable is the number of edges that were added during the
     * computation*/
    int num_edges;
} ripserResults;

template <typename DistanceMatrix>
class ripser
{
    const DistanceMatrix dist;
    index_t n, dim_max;
    const value_t threshold;
    const float ratio;
    const coefficient_t modulus;
    const binomial_coeff_table binomial_coeff;
    const std::vector<coefficient_t> multiplicative_inverse;
    mutable std::vector<diameter_entry_t> cofacet_entries;
    // If this flag is off, don't extract the representative cocycles to save
    // time
    const int do_cocycles;

    struct entry_hash {
        std::size_t operator()(const entry_t& e) const
        {
#if defined(USE_ROBINHOOD_HASHMAP)
            return robin_hood::hash<index_t>()(::get_index(e));
#else
            return std::hash<index_t>()(::get_index(e));
#endif
        }
    };

    struct equal_index {
        bool operator()(const entry_t& e, const entry_t& f) const
        {
            return ::get_index(e) == ::get_index(f);
        }
    };

    typedef hash_map<entry_t, size_t, entry_hash, equal_index> entry_hash_map;

public:
    mutable std::vector<std::vector<value_t>> births_and_deaths_by_dim;
    mutable std::vector<std::vector<std::vector<int>>> cocycles_by_dim;

    ripser(DistanceMatrix&& _dist, index_t _dim_max, value_t _threshold,
           float _ratio, coefficient_t _modulus, int _do_cocycles)
        : dist(std::move(_dist)), n(dist.size()), dim_max(_dim_max),
          threshold(_threshold), ratio(_ratio), modulus(_modulus),
          binomial_coeff(n, dim_max + 2),
          multiplicative_inverse(multiplicative_inverse_vector(_modulus)),
          do_cocycles(_do_cocycles)
    {
    }

    void copy_results(ripserResults& res)
    {
        res.births_and_deaths_by_dim = births_and_deaths_by_dim;
        res.cocycles_by_dim = cocycles_by_dim;
    }

    index_t get_max_vertex(const index_t idx, const index_t k,
                           const index_t n) const
    {
        auto top = n;
        auto bottom = k - 1;
        if (!(binomial_coeff(top, k) <= idx)) {
            index_t count = top - bottom;
            index_t step;
            index_t mid;
            while (count > 0) {
                step = count >> 1;
                mid = top - step;
                if (!(binomial_coeff(mid, k) <= idx)) {
                    top = mid - 1;
                    count -= step + 1;
                } else
                    count = step;
            }
        }
        return top;
    }

    index_t get_edge_index(const index_t i, const index_t j) const
    {
        return binomial_coeff(i, 2) + j;
    }

    template <typename OutputIterator>
    OutputIterator get_simplex_vertices(index_t idx, const index_t dim,
                                        index_t n, OutputIterator out) const
    {
        --n;
        for (index_t k = dim + 1; k > 0; --k) {
            n = get_max_vertex(idx, k, n);
            *out++ = n;
            idx -= binomial_coeff(n, k);
        }
        return out;
    }

    class simplex_coboundary_enumerator;

    void
    assemble_columns_to_reduce(std::vector<diameter_index_t>& simplices,
                               std::vector<diameter_index_t>& columns_to_reduce,
                               entry_hash_map& pivot_column_index, index_t dim)
    {
#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << "assembling columns" << std::flush;
        std::chrono::steady_clock::time_point next =
            std::chrono::steady_clock::now() + time_step;
#endif
        --dim;
        columns_to_reduce.clear();
        std::vector<diameter_index_t> next_simplices;

        for (diameter_index_t& simplex : simplices) {
            simplex_coboundary_enumerator cofacets(diameter_entry_t(simplex, 1),
                                                   dim, *this);

            while (cofacets.has_next(false)) {
#ifdef INDICATE_PROGRESS
                if (std::chrono::steady_clock::now() > next) {
                    std::cerr
                        << clear_line << "assembling " << next_simplices.size()
                        << " columns (processing "
                        << std::distance(&simplices[0], &simplex) << "/"
                        << simplices.size() << " simplices)" << std::flush;
                    next = std::chrono::steady_clock::now() + time_step;
                }
#endif
                auto cofacet = cofacets.next();
                if (get_diameter(cofacet) <= threshold) {
                    if (dim != dim_max)
                        next_simplices.push_back(
                            {get_diameter(cofacet), get_index(cofacet)});

                    if (pivot_column_index.find(get_entry(cofacet)) ==
                        pivot_column_index.end())
                        columns_to_reduce.push_back(
                            {get_diameter(cofacet), get_index(cofacet)});
                }
            }
        }

        simplices.swap(next_simplices);

#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << "sorting " << columns_to_reduce.size()
                  << " columns" << std::flush;
#endif

        std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
                  greater_diameter_or_smaller_index<diameter_index_t>());

#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << std::flush;
#endif
    }

    value_t get_vertex_birth(index_t i);
    void compute_dim_0_pairs(std::vector<diameter_index_t>& edges,
                             std::vector<diameter_index_t>& columns_to_reduce)
    {
        // TODO: Get correct birth times if the edges are negative (required for
        // lower star)
        union_find dset(n);
        for (index_t i = 0; i < n; i++) {
            dset.set_birth(i, get_vertex_birth(i));
        }

        edges = get_edges();
        std::sort(edges.rbegin(), edges.rend(),
                  greater_diameter_or_smaller_index<diameter_index_t>());
        std::vector<index_t> vertices_of_edge(2);
        for (auto e : edges) {
            get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.rbegin());
            index_t u = dset.find(vertices_of_edge[0]),
                    v = dset.find(vertices_of_edge[1]);

            if (u != v) {
                // Elder rule; youngest class (max birth time of u and v)
                // dies first
                value_t birth =
                    std::max(dset.get_birth(u), dset.get_birth(v));
                value_t death = get_diameter(e);
                if (death > birth) {
                    births_and_deaths_by_dim[0].push_back(birth);
                    births_and_deaths_by_dim[0].push_back(death);
                }
                dset.link(u, v);
            } else
                columns_to_reduce.push_back(e);
        }
        std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

        for (index_t i = 0; i < n; ++i)
            if (dset.find(i) == i) {
                births_and_deaths_by_dim[0].push_back(dset.get_birth(i));
                births_and_deaths_by_dim[0].push_back(
                    std::numeric_limits<value_t>::infinity());
            }
    }

    template <typename Column>
    diameter_entry_t pop_pivot(Column& column)
    {
        diameter_entry_t pivot(-1);
#ifdef USE_COEFFICIENTS
        while (!column.empty()) {
            if (get_coefficient(pivot) == 0)
                pivot = column.top();
            else if (get_index(column.top()) != get_index(pivot))
                return pivot;
            else
                set_coefficient(pivot,
                                get_modulo((get_coefficient(pivot) +
                                            get_coefficient(column.top())),
                                           modulus));
            column.pop();
        }
        return (get_coefficient(pivot) == 0) ? -1 : pivot;
#else
        while (!column.empty()) {
            pivot = column.top();
            column.pop();
            if (column.empty() || get_index(column.top()) != get_index(pivot))
                return pivot;
            column.pop();
        }
        return -1;
#endif
    }

    template <typename Column>
    diameter_entry_t get_pivot(Column& column)
    {
        diameter_entry_t result = pop_pivot(column);
        if (get_index(result) != -1)
            column.push(result);
        return result;
    }

    template <typename Column>
    diameter_entry_t init_coboundary_and_get_pivot(
        const diameter_entry_t simplex, Column& working_coboundary,
        const index_t& dim, entry_hash_map& pivot_column_index)
    {
        bool check_for_emergent_pair = true;
        cofacet_entries.clear();
        simplex_coboundary_enumerator cofacets(simplex, dim, *this);
        while (cofacets.has_next()) {
            diameter_entry_t cofacet = cofacets.next();
            if (get_diameter(cofacet) <= threshold) {
                cofacet_entries.push_back(cofacet);
                if (check_for_emergent_pair &&
                    (get_diameter(simplex) == get_diameter(cofacet))) {
                    if (pivot_column_index.find(get_entry(cofacet)) ==
                        pivot_column_index.end())
                        return cofacet;
                    check_for_emergent_pair = false;
                }
            }
        }
        for (auto cofacet : cofacet_entries)
            working_coboundary.push(cofacet);
        return get_pivot(working_coboundary);
    }

    template <typename Column>
    void add_simplex_coboundary(const diameter_entry_t simplex,
                                const index_t& dim,
                                Column& working_reduction_column,
                                Column& working_coboundary)
    {
        working_reduction_column.push(simplex);
        simplex_coboundary_enumerator cofacets(simplex, dim, *this);
        while (cofacets.has_next()) {
            diameter_entry_t cofacet = cofacets.next();
            if (get_diameter(cofacet) <= threshold)
                working_coboundary.push(cofacet);
        }
    }

    template <typename Column>
    void
    add_coboundary(compressed_sparse_matrix<diameter_entry_t>& reduction_matrix,
                   const std::vector<diameter_index_t>& columns_to_reduce,
                   const size_t index_column_to_add, const coefficient_t factor,
                   const size_t& dim, Column& working_reduction_column,
                   Column& working_coboundary)
    {
        diameter_entry_t column_to_add(columns_to_reduce[index_column_to_add],
                                       factor);
        add_simplex_coboundary(column_to_add, dim, working_reduction_column,
                               working_coboundary);

        for (diameter_entry_t simplex :
             reduction_matrix.subrange(index_column_to_add)) {
            set_coefficient(simplex,
                            get_coefficient(simplex) * factor % modulus);
            add_simplex_coboundary(simplex, dim, working_reduction_column,
                                   working_coboundary);
        }
    }

    using working_t = std::priority_queue<
        diameter_entry_t, std::vector<diameter_entry_t>,
        greater_diameter_or_smaller_index<diameter_entry_t>>;

    diameter_entry_t cocycle_e;
    std::vector<index_t> cocycle_simplex;
    std::vector<int> thiscocycle;
    inline void compute_cocycles(working_t cocycle, index_t dim)
    {
        thiscocycle.clear();
        while (get_index(cocycle_e = get_pivot(cocycle)) != -1) {
            cocycle_simplex.clear();
            get_simplex_vertices(get_index(cocycle_e), dim, n,
                                 std::back_inserter(cocycle_simplex));
            for (size_t k = 0; k < cocycle_simplex.size(); k++) {
                thiscocycle.push_back((int) cocycle_simplex[k]);
            }
            thiscocycle.push_back(
                normalize(get_coefficient(cocycle_e), modulus));
            cocycle.pop();
        }
        cocycles_by_dim[dim].push_back(thiscocycle);
    }

    void compute_pairs(std::vector<diameter_index_t>& columns_to_reduce,
                       entry_hash_map& pivot_column_index, index_t dim)
    {
        compressed_sparse_matrix<diameter_entry_t> reduction_matrix;
        size_t index_column_to_add;

#ifdef INDICATE_PROGRESS
        std::chrono::steady_clock::time_point next =
            std::chrono::steady_clock::now() + time_step;
#endif

        for (size_t index_column_to_reduce = 0;
             index_column_to_reduce < columns_to_reduce.size();
             ++index_column_to_reduce) {
            diameter_entry_t column_to_reduce(
                columns_to_reduce[index_column_to_reduce], 1);
            value_t diameter = get_diameter(column_to_reduce);

            reduction_matrix.append_column();

            working_t working_reduction_column;
            working_t working_coboundary;

            working_reduction_column.push(column_to_reduce);

            diameter_entry_t pivot = init_coboundary_and_get_pivot(
                column_to_reduce, working_coboundary, dim, pivot_column_index);

            while (true) {
#ifdef INDICATE_PROGRESS
                if (std::chrono::steady_clock::now() > next) {
                    std::cerr << clear_line << "reducing column "
                              << index_column_to_reduce + 1 << "/"
                              << columns_to_reduce.size() << " (diameter "
                              << diameter << ")" << std::flush;
                    next = std::chrono::steady_clock::now() + time_step;
                }
#endif
                if (get_index(pivot) != -1) {
                    auto pair = pivot_column_index.find(get_entry(pivot));
                    if (pair != pivot_column_index.end()) {
                        entry_t other_pivot = pair->first;
                        index_column_to_add = pair->second;
                        coefficient_t factor =
                            modulus -
                            get_modulo(
                                get_coefficient(pivot) *
                                    multiplicative_inverse[get_coefficient(
                                        other_pivot)],
                                modulus);

                        add_coboundary(reduction_matrix, columns_to_reduce,
                                       index_column_to_add, factor, dim,
                                       working_reduction_column,
                                       working_coboundary);

                        pivot = get_pivot(working_coboundary);
                    } else {
                        value_t death = get_diameter(pivot);
                        if (death > diameter * ratio) {
                            births_and_deaths_by_dim[dim].push_back(diameter);
                            births_and_deaths_by_dim[dim].push_back(death);
                            if (do_cocycles) {
                                // Representative cocycle
                                compute_cocycles(working_reduction_column, dim);
                            }
                        }

                        pivot_column_index.insert(
                            {get_entry(pivot), index_column_to_reduce});

                        pop_pivot(working_reduction_column);
                        while (true) {
                            diameter_entry_t e =
                                pop_pivot(working_reduction_column);

                            if (get_index(e) == -1)
                                break;
                            assert(get_coefficient(e) > 0);
                            reduction_matrix.push_back(e);
                        }
                        break;
                    }
                } else {
                    births_and_deaths_by_dim[dim].push_back(diameter);
                    births_and_deaths_by_dim[dim].push_back(
                        std::numeric_limits<value_t>::infinity());

                    if (do_cocycles) {
                        // Representative cocycle
                        compute_cocycles(working_reduction_column, dim);
                    }
                    break;
                }
            }
        }
#ifdef INDICATE_PROGRESS
        std::cerr << clear_line << std::flush;
#endif
    }

    std::vector<diameter_index_t> get_edges();
    void compute_barcodes()
    {
        std::vector<diameter_index_t> simplices, columns_to_reduce;

        /* prevent cases where dim_max < 0 */
        if (dim_max < 0)
            dim_max = 0;

        births_and_deaths_by_dim.resize(dim_max + 1);
        cocycles_by_dim.resize(dim_max + 1);

        compute_dim_0_pairs(simplices, columns_to_reduce);

        for (index_t dim = 1; dim <= dim_max; ++dim) {
            entry_hash_map pivot_column_index;
            pivot_column_index.reserve(columns_to_reduce.size());

            compute_pairs(columns_to_reduce, pivot_column_index, dim);

            if (dim < dim_max)
                assemble_columns_to_reduce(simplices, columns_to_reduce,
                                           pivot_column_index, dim + 1);
        }
    }
};

template <>
value_t ripser<compressed_lower_distance_matrix>::get_vertex_birth(index_t i)
{
    // TODO: Dummy for now; nonzero vertex births are only done through
    // sparse matrices at the moment
    return 0.0;
}

template <>
value_t ripser<sparse_distance_matrix>::get_vertex_birth(index_t i)
{
    return dist.vertex_births[i];
}

template <>
class ripser<compressed_lower_distance_matrix>::simplex_coboundary_enumerator
{
private:
    index_t idx_below, idx_above, v, k;
    std::vector<index_t> vertices;
    const diameter_entry_t simplex;
    const coefficient_t modulus;
    const compressed_lower_distance_matrix& dist;
    const binomial_coeff_table& binomial_coeff;

public:
    simplex_coboundary_enumerator(
        const diameter_entry_t _simplex, index_t _dim,
        const ripser<compressed_lower_distance_matrix>& parent)
        : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1),
          k(_dim + 1), vertices(_dim + 1), simplex(_simplex),
          modulus(parent.modulus), dist(parent.dist),
          binomial_coeff(parent.binomial_coeff)
    {
        parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n,
                                    vertices.begin());
    }

    bool has_next(bool all_cofacets = true)
    {
        return (v >= k && (all_cofacets || binomial_coeff(v, k) > idx_below));
    }

    diameter_entry_t next()
    {
        while ((binomial_coeff(v, k) <= idx_below)) {
            idx_below -= binomial_coeff(v, k);
            idx_above += binomial_coeff(v, k + 1);
            --v;
            --k;
            assert(k != -1);
        }
        value_t cofacet_diameter = get_diameter(simplex);
        for (index_t w : vertices)
            cofacet_diameter = std::max(cofacet_diameter, dist(v, w));
        index_t cofacet_index =
            idx_above + binomial_coeff(v--, k + 1) + idx_below;
        coefficient_t cofacet_coefficient =
            (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;
        return diameter_entry_t(cofacet_diameter, cofacet_index,
                                cofacet_coefficient);
    }
};

template <>
class ripser<sparse_distance_matrix>::simplex_coboundary_enumerator
{
    index_t idx_below, idx_above, k;
    std::vector<index_t> vertices;
    const diameter_entry_t simplex;
    const coefficient_t modulus;
    const sparse_distance_matrix& dist;
    const binomial_coeff_table& binomial_coeff;
    std::vector<std::vector<index_diameter_t>::const_reverse_iterator>&
        neighbor_it;
    std::vector<std::vector<index_diameter_t>::const_reverse_iterator>&
        neighbor_end;
    index_diameter_t neighbor;

public:
    simplex_coboundary_enumerator(const diameter_entry_t _simplex,
                                  const index_t _dim,
                                  const ripser<sparse_distance_matrix>& parent)
        : idx_below(get_index(_simplex)), idx_above(0), k(_dim + 1),
          vertices(_dim + 1), simplex(_simplex), modulus(parent.modulus),
          dist(parent.dist), binomial_coeff(parent.binomial_coeff),
          neighbor_it(dist.neighbor_it), neighbor_end(dist.neighbor_end)
    {
        neighbor_it.clear();
        neighbor_end.clear();

        parent.get_simplex_vertices(idx_below, _dim, parent.n,
                                    vertices.rbegin());

        for (auto v : vertices) {
            neighbor_it.push_back(dist.neighbors[v].rbegin());
            neighbor_end.push_back(dist.neighbors[v].rend());
        }
    }

    bool has_next(bool all_cofacets = true)
    {
        for (auto &it0 = neighbor_it[0], &end0 = neighbor_end[0]; it0 != end0;
             ++it0) {
            neighbor = *it0;
            for (size_t idx = 1; idx < neighbor_it.size(); ++idx) {
                auto &it = neighbor_it[idx], end = neighbor_end[idx];
                while (get_index(*it) > get_index(neighbor))
                    if (++it == end)
                        return false;
                if (get_index(*it) != get_index(neighbor))
                    goto continue_outer;
                else
                    neighbor = std::max(neighbor, *it);
            }
            while (k > 0 && vertices[k - 1] > get_index(neighbor)) {
                if (!all_cofacets)
                    return false;
                idx_below -= binomial_coeff(vertices[k - 1], k);
                idx_above += binomial_coeff(vertices[k - 1], k + 1);
                --k;
            }
            return true;
        continue_outer:;
        }
        return false;
    }

    diameter_entry_t next()
    {
        ++neighbor_it[0];
        value_t cofacet_diameter =
            std::max(get_diameter(simplex), get_diameter(neighbor));
        index_t cofacet_index =
            idx_above + binomial_coeff(get_index(neighbor), k + 1) + idx_below;
        coefficient_t cofacet_coefficient =
            (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;
        return diameter_entry_t(cofacet_diameter, cofacet_index,
                                cofacet_coefficient);
    }
};

template <>
std::vector<diameter_index_t>
ripser<compressed_lower_distance_matrix>::get_edges()
{
    std::vector<diameter_index_t> edges;
    std::vector<index_t> vertices(2);
    for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
        get_simplex_vertices(index, 1, dist.size(), vertices.rbegin());
        value_t length = dist(vertices[0], vertices[1]);
        if (length <= threshold)
            edges.push_back({length, index});
    }
    return edges;
}

template <>
std::vector<diameter_index_t> ripser<sparse_distance_matrix>::get_edges()
{
    std::vector<diameter_index_t> edges;
    for (index_t i = 0; i < n; ++i)
        for (auto n : dist.neighbors[i]) {
            index_t j = get_index(n);
            if (i > j)
                edges.push_back({get_diameter(n), get_edge_index(i, j)});
        }
    return edges;
}

ripserResults rips_dm(float* D, int N, int modulus, int dim_max,
                      float threshold, int do_cocycles)
{
    // Setup distance matrix and figure out threshold
    std::vector<value_t> distances(D, D + N);
    compressed_lower_distance_matrix dist = compressed_lower_distance_matrix(
        compressed_upper_distance_matrix(std::move(distances)));

    // TODO: This seems like a dummy parameter at the moment
    float ratio = 1.0;

    value_t min = std::numeric_limits<value_t>::infinity(),
            max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
    int num_edges = 0;

    /* Use enclosing radius when users does not set threshold or
     * when users uses infinity as a threshold
     */
    if (threshold == std::numeric_limits<value_t>::max() ||
        threshold == std::numeric_limits<value_t>::infinity()) {
        value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
        for (size_t i = 0; i < dist.size(); ++i) {
            value_t r_i = -std::numeric_limits<value_t>::infinity();
            for (size_t j = 0; j < dist.size(); ++j)
                r_i = std::max(r_i, dist(i, j));
            enclosing_radius = std::min(enclosing_radius, r_i);
        }
        threshold = enclosing_radius;
    }

    for (auto d : dist.distances) {
        min = std::min(min, d);
        max = std::max(max, d);
        max_finite = d != std::numeric_limits<value_t>::infinity()
                         ? std::max(max, d)
                         : max_finite;
        if (d <= threshold)
            ++num_edges;
    }

    ripserResults res;
    if (threshold >= max) {
        ripser<compressed_lower_distance_matrix> r(
            std::move(dist), dim_max, threshold, ratio, modulus, do_cocycles);
        r.compute_barcodes();
        r.copy_results(res);
    } else {
        ripser<sparse_distance_matrix> r(
            sparse_distance_matrix(std::move(dist), threshold), dim_max,
            threshold, ratio, modulus, do_cocycles);
        r.compute_barcodes();
        r.copy_results(res);
    }
    res.num_edges = num_edges;
    return res;
}

ripserResults rips_dm_sparse(int* I, int* J, float* V, int NEdges, int N,
                             int modulus, int dim_max, float threshold,
                             int do_cocycles)
{
    // TODO: This seems like a dummy parameter at the moment
    float ratio = 1.0;
    // Setup distance matrix and figure out threshold
    ripser<sparse_distance_matrix> r(
        sparse_distance_matrix(I, J, V, NEdges, N, threshold), dim_max,
        threshold, ratio, modulus, do_cocycles);
    r.compute_barcodes();
    // Report the number of edges that were added
    int num_edges = 0;
    for (int idx = 0; idx < NEdges; idx++) {
        if (I[idx] < J[idx] && V[idx] <= threshold) {
            num_edges++;
        }
    }
    ripserResults res;
    r.copy_results(res);
    res.num_edges = num_edges;
    return res;
}

#ifdef LIBRIPSER
int unpack_results(int** n_intervals, value_t** births_and_deaths,
                   int** cocycle_length, int** cocycles, ripserResults res,
                   int do_cocycles)
{
    int n_dims = res.births_and_deaths_by_dim.size();
    *n_intervals = (int*) malloc(n_dims * sizeof(int));
    int n_intervals_total = 0;

    for (int d = 0; d < n_dims; d++) {
        int n_int_d = res.births_and_deaths_by_dim[d].size() / 2;
        (*n_intervals)[d] = n_int_d;
        n_intervals_total += n_int_d;
    }
    *births_and_deaths =
        (value_t*) malloc(n_intervals_total * 2 * sizeof(value_t));
    *cocycle_length = (int*) calloc(n_intervals_total, sizeof(int));

    int cocycle_length_total = 0;
    int idx = 0;
    for (int d = 0; d < n_dims; d++) {
        std::copy(res.births_and_deaths_by_dim[d].begin(),
                  res.births_and_deaths_by_dim[d].end(),
                  &(*births_and_deaths)[2 * idx]);

        if (do_cocycles && !res.cocycles_by_dim[d].empty()) {
            for (int i = 0; i < (*n_intervals)[d]; i++) {
                int cc_length = res.cocycles_by_dim[d][i].size();
                (*cocycle_length)[idx] = cc_length;
                cocycle_length_total += cc_length;
                idx++;
            }
        } else {
            idx += (*n_intervals)[d];
        }
    }

    if (do_cocycles && cocycle_length_total > 0) {
        *cocycles = (int*) calloc(cocycle_length_total, sizeof(int));

        int pos = 0;
        for (int d = 0; d < n_dims; d++) {
            if (!res.cocycles_by_dim[d].empty()) {
                for (int i = 0; i < (*n_intervals)[d]; i++) {
                    int cc_length = res.cocycles_by_dim[d][i].size();
                    std::copy(res.cocycles_by_dim[d][i].begin(),
                              res.cocycles_by_dim[d][i].end(),
                              &(*cocycles)[pos]);
                    pos += cc_length;
                }
            }
        }
    }
    return res.num_edges;
}
extern "C" {
#include "ripser.h"

/*
  C interface to Ripser.

  Results are passed through output arguments. The arrays are allocated in this
  function and have to be freed manually by the caller.

  Output arguments:
  * n_intervals: number of intervals per dimension. (length = dim_max + 1)
  * births_and_deaths: births and deaths of all dimension in a flat array.
  (length = 2 * sum(n_intervals))
  * cocycle_length: lengths of individual cocycles. (length = sum(n_intervals))
  * cocycles: cocycles stored in a flat array. (length = sum(cocycle_length))
  Input arguments:
  * D: lower triangle of the distance matrix in a flat array.
  * N: length of D.
  * modulus: Compute homology with coefficients in the prime field Z/pZ. p must
  be a prime number.
  * dim_max: Compute persistent homology up to this dimension
  * threshold: Compute Rips complexes up to this diameter
  * do_cocycles: If nonzero, calculate cocycles and write them to cocycle_length
  and cocycles.

  Returns number of edges.
*/
int c_rips_dm(int** n_intervals, value_t** births_and_deaths,
              int** cocycle_length, int** cocycles, value_t* D, int N,
              int modulus, int dim_max, value_t threshold, int do_cocycles)
{
    ripserResults res = rips_dm(D, N, modulus, dim_max, threshold, do_cocycles);
    return unpack_results(n_intervals, births_and_deaths, cocycle_length,
                          cocycles, res, do_cocycles);
}

/*
  Same as c_rips_dm, but takes sparse matrix as input.

  Arguments:
  * I, J, V: indices and values.
  * NEdges: total number of indices and values.
  * N: size of sparse matrix.
*/
int c_rips_dm_sparse(int** n_intervals, value_t** births_and_deaths,
                     int** cocycle_length, int** cocycles, int* I, int* J,
                     float* V, int NEdges, int N, int modulus, int dim_max,
                     float threshold, int do_cocycles)
{
    ripserResults res = rips_dm_sparse(I, J, V, NEdges, N, modulus, dim_max,
                                       threshold, do_cocycles);
    return unpack_results(n_intervals, births_and_deaths, cocycle_length,
                          cocycles, res, do_cocycles);
}
}

#endif

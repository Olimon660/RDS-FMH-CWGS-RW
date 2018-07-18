/*
 * ./dream dataset.txt package1.txt output1.txt
 */

#include <set>
#include <map>
#include <memory>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "phone.hpp"

using namespace DreamLab;

typedef unsigned Index;
typedef std::string Line;
typedef std::string FileName;
typedef std::set<std::size_t> Missing;

struct Package
{
    std::set<std::size_t> NA;
    
    // Number of permutations
    unsigned perm;
    
    std::string mode;
    
    Seed s1, s2, s3;
    
    // Random package ID
    std::string id;
    
    // Matrix for phenotype
    std::shared_ptr<Matrix> Y;
    
    // Matrix for covariates and genotype
    std::shared_ptr<Matrix> X;
};

struct Data
{
    // Random data ID
    std::string id;

    std::set<std::size_t> NA;
    std::vector<long double> G;
};

static std::shared_ptr<Matrix> filter(const Matrix &src, const Missing &m)
{
    auto dst = std::shared_ptr<Matrix>(new Matrix(src.rows() - m.size(), src.cols()));
    auto i_ = 0;
    
    for (auto i = 0; i < src.rows(); i++)
    {
        if (!m.count(i))
        {
            for (auto j = 0; j < src.cols(); j++)
            {
                (*dst)(i_, j) = src(i, j);
            }

            i_++;
        }
    }
    
    return dst;
}

template<typename T> void split(const std::string &s, char delim, T result)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

static Package readP(const FileName &file)
{
    std::ifstream f(file);
    
    Line l;
    std::getline(f, l);

    std::vector<std::string> toks;
	split(l, '#', std::back_inserter(toks));
        
    Package p;
    p.id = toks[0];
    p.mode = toks[4];
    p.perm = stoi(toks[1]);

    std::vector<std::string> seeds;
	split(toks[5], ',', std::back_inserter(seeds));
    p.s1 = stoi(seeds[0]); p.s2 = stoi(seeds[1]); p.s3 = stoi(seeds[2]);
    assert(p.s1 && p.s2 && p.s3);

    std::cout << "Number of permutations: " << p.perm << std::endl;
        
    const auto ns = stoi(toks[2]);
    const auto nc = stoi(toks[3]);

    p.Y = std::shared_ptr<Matrix>(new Matrix(ns, 1));
    p.X = std::shared_ptr<Matrix>(new Matrix(ns, 1 + nc + 1));
    p.X->setOnes();
        
    // Which column does the variables start?
    const auto start = 6;

    for (auto j = start; j < toks.size(); j++)
    {
        std::vector<std::string> vs;
		split(toks[j], ',', std::back_inserter(vs));
        
        if (j == start) // Phenotype?
        {
            for (auto k = 0; k < vs.size(); k++)
            {
                if (vs[k] == "nan") { p.NA.insert(k); (*p.Y)(k, 0) = NAN; }
                else                {                 (*p.Y)(k, 0) = stold(vs[k]); }
            }
        }
        else
        {
            for (auto k = 0; k < vs.size(); k++)
            {
                if (vs[k] == "nan") { p.NA.insert(k); (*p.X)(k, 0) = NAN; }
                else                {                 (*p.X)(k, j-start+1) = stold(vs[k]); }
            }
        }
    }

    return p;
}

static Data readG(const std::vector<std::string> &data)
{
    std::vector<std::string> g;
	split(data[1], ',', std::back_inserter(g));

    Data x;
    
    for (auto i = 0; i < g.size(); i++)
    {
        if (g[i] == "3") { x.NA.insert(i); x.G.push_back(NAN); }
        else             { x.G.push_back(stof(g[i])); }
    }
    
    return x;
}

#ifdef DEBUG
static void writeM(const FileName &file, const Matrix &m)
{
    std::ofstream f(file);
    if (f.is_open())
    {
        for (auto i = 0; i < m.rows(); i++)
        {
            for (auto j = 0; j < m.cols(); j++)
            {
                if (j) { f << ","; } f << m(i,j);
            }
            
            f << "\n";
        }
    }
    f.close();
}
#endif

static void compute(const Data &data, const Package &pack, std::stringstream &ss)
{
    std::cout << "Number of samples filtered from data: " << data.NA.size() << std::endl;

    Missing NA;
    NA.insert(pack.NA.begin(), pack.NA.end());
    NA.insert(data.NA.begin(), data.NA.end());

    auto X_ = filter(*(pack.X), NA);
    auto Y_ = filter(*(pack.Y), NA);

    auto j = 0;
    for (auto i = 0; i < data.G.size(); i++)
    {
        if (!NA.count(i))
        {
            assert(data.G[i] != 3);
            (*X_)(j++, 1) = data.G[i];
        }
    }
    
    const auto X = *X_;
    const auto Y = *Y_;
  
    const auto B = pack.perm;
    const auto mode = pack.mode;
    
#ifdef DEBUG
    writeM("X.txt", X);
    writeM("Y.txt", Y);
#endif

    // Results without permutation
    GStats gs = (mode == "linear") ? _gt_statistic<Linear>(Y, X) : (mode == "logistic") ? _gt_statistic<Logistic>(Y, X) : _gt_statistic<Firth>(Y, X);
    
    auto toString = [](const GStats &x)
    {
        #define STR(x) std::to_string(x)
        return "(" + STR(std::get<0>(x)) + "," + STR(std::get<1>(x)) + "," + STR(std::get<2>(x)) + ")";
    };
    
    const auto s1 = pack.s1;
    const auto s2 = pack.s2;
    const auto s3 = pack.s3;
    
    // Results after permuation
    const auto bs = (mode == "linear") ? _gt_statistic_null<Linear>(Y, X, B, s1, s2, s3) : (mode == "logistic") ? _gt_statistic_null<Logistic>(Y, X, B, s1, s2, s3) : _gt_statistic_null<Firth>(Y, X, B, s1, s2, s3);
    
    ss << data.id << "#" << pack.id << "#" << toString(gs) << "#";
    for (auto i = 0; i < bs.size(); i++) { if (i != 0) { ss << ","; } ss << toString(bs[i]); }
    ss << "\n";
}

void dreamLab(const FileName &i1, const FileName &i2, const FileName &o)
{
    const auto P = readP(i2);
    std::cout << "Number of samples filtered from package: " << P.NA.size() << std::endl;

    Line l;
    std::stringstream ss;

    unsigned i = 0;
    std::ifstream f1(i1);
    
    while (std::getline(f1, l))
    {
        std::cout << i++ << std::endl;

        std::vector<std::string> data;
		split(l, '#', std::back_inserter(data));

        compute(readG(data), P, ss);
    }

    std::ofstream w;
    w.open(o);
    w << ss.str();
    w.close();
}

#ifndef UNIT_TEST
int main(int argc, char** argv)
{
    dreamLab(argv[1], argv[2], argv[3]);
    return 0;
}
#endif

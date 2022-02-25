
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <ranges>
#include <filesystem>

#include "common.hpp"
#include "primitives.hpp"
#include "aabb.hpp"

#include "serializer.hpp"
#include "print_progress.hpp"

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>


namespace mirrors {

enum class source_type {
    directional,
    point,
};

using mirror_t = std::vector<line_t>;
struct mirrors_t {
    mirror_t m1,m2;
    point_t centre;
    f_t D;
    aabb<f_t> bb;
    source_type src_type;
    point_t src;
    f_t point_src_max_alpha{};
    std::vector<std::pair<float,float>> R_of_beta;
};

inline mirror_t construct_mirror(const mirror_t &m, const point_t &flip) {
    mirror_t ret{};
    ret.reserve(m.size()*2);
    std::ranges::reverse_view rv{m};
    for (const auto &s : rv)
        ret.emplace_back(s.a*flip,s.b*flip,s.n1*flip,s.n2*flip,s.mirror);
    ret.insert(ret.end(),m.begin(),m.end());
    return ret;
}

inline aabb<f_t> mirror_aabb(const mirror_t &m) {
    if (!m.size()) return {};
    aabb<f_t> bb = m.back().AABB();
    for (const auto &l : m) 
        bb += l.AABB();
    bb.b += point_t{1e-8};
    return bb;
}

template <typename Ifun, typename Sfun, typename ProgFun>
void compute_mirrors_design1_segment(const f_t f, const f_t L, const f_t beam_in_angle, bool invert, Ifun &&I, Sfun &&S, ProgFun &&progf, f_t dY, const f_t beta_max, mirror_t &m1, mirror_t &m2,
    std::vector<std::pair<float,float>> *R_of_beta) {
    f_t X{-L}, Y{}, 
        Xt{f}, Yt{}, 
        dX_dY{}, dXt_dYt{},
        b{}, r{f}, br{};
    point_t n, nt{ -1,0,0 };

    const auto& Kin = -point_t{ cos(beam_in_angle),sin(beam_in_angle),0 };
    n = point_t{ cos(beam_in_angle/2),sin(beam_in_angle/2),0 };
    dY *= cos(beam_in_angle/2);
    dX_dY = -n.y/n.x;

    for (;glm::abs(b)<=beta_max;) {
        const auto p  = glm::abs(b)/beta_max;

        const auto Ynew  = Y + dY;
        const auto dX    = dX_dY*(Ynew-Y);
        const auto Xnew  = X + dX;
        const auto brnew = glm::length(glm::cross(Kin, point_t{ Xnew+L,Ynew,0 }));
        const auto db    = (invert?1:-1) * I(br)/S(b) * (brnew-br) * glm::sign(dY);
        const auto bnew  = b + db;

        const auto step  = point_t{ nt.y,-nt.x,0 };
        /* if (step.y<=0) { */
        /*     std::cerr << "Maximal possible beta (=" << std::to_string(glm::abs(b)*180/glm::pi<f_t>()) << ") reached" << std::endl; */
        /*     break; */
        /* } */
        const auto phat  = point_t{ glm::cos(b),glm::sin(b),0 };
        const auto sina  = glm::length(glm::cross(phat,step));
        const auto pnew  = point_t{ Xt,Yt,0 } + r*glm::sin(db)/sina * step;

        const auto Ytnew = pnew.y;
        const auto Xtnew = pnew.x;
        /* const auto Ytnew = r*glm::sin(bnew); */
        /* const auto Xtnew = Xt + dXt_dYt*(Ytnew-Yt); */

        const f_t  rnew  = glm::sqrt(sqr(Xtnew) + sqr(Ytnew));
        const f_t  R     = glm::sqrt(sqr(Xtnew-Xnew) + sqr(Ytnew-Ynew));
	if (R_of_beta)
		R_of_beta->emplace_back(glm::abs(bnew),R+rnew);

        const auto Kint  = f_t(1)/R * point_t{ Xtnew-Xnew, Ytnew-Ynew, 0 };
        const auto Kout  = -point_t{ glm::cos(bnew),glm::sin(bnew),0 };
        const auto Kn1   = Kint-Kin;
        const auto Kn2   = Kout-Kint;

        const auto dX_dY_new   = -Kn1.y/Kn1.x;
        const auto dXt_dYt_new = -Kn2.y/Kn2.x;

        const auto nnew  = glm::normalize(Kn1);
        const auto ntnew = glm::normalize(Kn2);
        m1.emplace_back(point_t{ X,  Y,  0 }, point_t{ Xnew,  Ynew,  0 }, n,  nnew,  1);
        m2.emplace_back(point_t{ Xt, Yt, 0 }, point_t{ Xtnew, Ytnew, 0 }, nt, ntnew, 2);

        X  = Xnew;
        Y  = Ynew;
        Xt = Xtnew;
        Yt = Ytnew;
        r  = rnew;
        b  = bnew;
        br = brnew;
        dX_dY   = dX_dY_new;
        dXt_dYt = dXt_dYt_new;
        n  = nnew;
        nt = ntnew;

        progf(p);
    }
}
template <typename Ifun, typename Sfun>
mirrors_t compute_mirrors_design1(const f_t f, const f_t L, const f_t beam_in_angle, bool invert, Ifun &&I, Sfun &&S, const f_t dY, const f_t beta_max) {
    mirror_t m1l{}, m2l{};
    mirror_t m1r{}, m2r{};

    std::vector<std::pair<float,float>> R_of_beta;
    
    compute_mirrors_design1_segment(f,L,beam_in_angle,invert,I,S,[](auto p) { print_progress(p/2);    },-dY,beta_max,m1l,m2l,&R_of_beta);
    compute_mirrors_design1_segment(f,L,beam_in_angle,invert,I,S,[](auto p) { print_progress(.5+p/2); }, dY,beta_max,m1r,m2r,nullptr);

    mirrors_t m;

    m.m1.reserve(m1l.size());
    std::ranges::reverse_view rv1{m1l}, rv2{m2l};
    for (const auto &s : rv1)
        m.m1.emplace_back(s);
    for (const auto &s : m1r)
        m.m1.emplace_back(s);
    for (const auto &s : rv2)
        m.m2.emplace_back(s);
    for (const auto &s : m2r)
        m.m2.emplace_back(s);
    
    m.bb = mirror_aabb(m.m1) + mirror_aabb(m.m2);
    m.D  = glm::length(m.m1.back().b-m.m1.front().a);
    m.centre  = point_t{ -L,0,0 };
    m.src_type = source_type::directional;
    m.src = -point_t{ cos(beam_in_angle),sin(beam_in_angle),0 };
    m.R_of_beta = std::move(R_of_beta);

    return m;
}

template <typename Ifun, typename Sfun>
mirrors_t compute_mirrors_design3(const f_t f1, const f_t f2, const f_t L, Ifun &&S1, Sfun &&S2, const f_t dY, const f_t beta_max) {
    f_t X{-f1-L}, Y{}, 
        Xt{f2}, Yt{}, 
        dX_dY{}, dXt_dYt{},
        a{}, b{}, 
        r1{f1}, r2{f2};
    point_t n{ 1,0,0 },
            nt{ -1,0,0 };

    mirror_t m1{}, m2{};
    for (;b<=beta_max;) {
        const auto p  = b/beta_max;

        const auto Ynew  = Y + dY;
        const auto dX    = dX_dY * dY;
        const auto Xnew  = X + dX;
        const f_t  r1new = glm::sqrt(sqr(Xnew+L)+sqr(Ynew));
        const auto anew  = glm::asin(Ynew/r1new);
        if (anew<=a) {
            std::cerr << "beta_max (m) too large for setup" << std::endl;
            exit(1);
            return {};
        }
        const auto da    = anew - a;
        const auto db    = S1(a)/S2(b) * da;
        const auto bnew  = b + db;
        const auto Ytnew = r2 * glm::sin(bnew);
        const auto Xtnew = Xt + dXt_dYt * (Ytnew - Yt);
        const f_t  r2new = glm::sqrt(sqr(Xtnew) + sqr(Ytnew));
        const f_t  R     = glm::sqrt(sqr(Xtnew - Xnew) + sqr(Ytnew - Ynew));

        const auto Kin   = point_t{ -glm::cos(anew),glm::sin(anew),0 };
        const auto Kint  = f_t(1)/R * point_t{ Xtnew-Xnew, Ytnew-Ynew, 0 };
        const auto Kout  = -point_t{ glm::cos(bnew),glm::sin(bnew),0 };
        const auto Kn1   = Kint-Kin;
        const auto Kn2   = Kout-Kint;

        const auto dX_dY_new   = -Kn1.y/Kn1.x;
        const auto dXt_dYt_new = -Kn2.y/Kn2.x;

        const auto nnew  = glm::normalize(Kn1);
        const auto ntnew = glm::normalize(Kn2);
        m1.emplace_back(point_t{X,Y,0},   point_t{Xnew,  Ynew,0},  n,  nnew,  1);
        m2.emplace_back(point_t{Xt,Yt,0}, point_t{Xtnew, Ytnew,0}, nt, ntnew, 2);

        X  = Xnew;
        Y  = Ynew;
        Xt = Xtnew;
        Yt = Ytnew;
        r1 = r1new;
        r2 = r2new;
        a  = anew;
        b  = bnew;
        dX_dY   = dX_dY_new;
        dXt_dYt = dXt_dYt_new;
        n  = nnew;
        nt = ntnew;
        
        print_progress(p-1e-10);
    }

    mirrors_t m;
    m.m1 = construct_mirror(m1,point_t{1,-1,1});
    m.m2 = construct_mirror(m2,point_t{1,-1,1});
    m.bb = mirror_aabb(m.m1) + mirror_aabb(m.m2);
    m.bb.b += point_t{1e-5,1e-5,1e-5};
    m.D  = 2*m1.back().b.y;
    m.centre  = point_t{ -L,0,0 };
    m.src_type = source_type::point;
    m.src = m.centre;
    m.point_src_max_alpha = a;

    return m;
}

inline mirrors_t compute_mirror_parabolic(const f_t f, const f_t dY, const f_t beta_max) {
    mirror_t ml{},mr{};

    std::vector<std::pair<float,float>> R_of_beta;
    
    const auto inv = point_t{1,-1,1};
    point_t p{ -f,0,0 }, n{ 1,0,0 };
    f_t b{};
    for (;b<=beta_max;) {
        const auto y = p.y+dY;
        const auto x = sqr(p.y)/(4*f)-f;
        const auto nnew = glm::normalize(point_t{ 1,-y/2/f,0 });

        mr.emplace_back(p,     point_t{x, y,0}, n,     nnew,     2);
        ml.emplace_back(p*inv, point_t{x,-y,0}, n*inv, nnew*inv, 2);
        
        p = point_t{ x,y,0 };
        b = glm::atan(y,-x);
        n = nnew;

	R_of_beta.emplace_back(b,glm::length(p));
        
        print_progress(b/beta_max);
    }

    mirrors_t m;

    m.m2.reserve(ml.size()*2);
    std::ranges::reverse_view rv{ml};
    for (const auto &s : rv)
        m.m2.emplace_back(s);
    for (const auto &s : mr)
        m.m2.emplace_back(s);
    
    m.bb = aabb{ point_t {-f,-p.y,0 }, point_t{ p.x,p.y,1e-5 } };
    m.D  = 2*p.y;
    m.centre  = point_t{ 0,0,0 };
    m.src_type = source_type::directional;
    m.src = -point_t{ 1,0,0 };
    m.R_of_beta = std::move(R_of_beta);

    return m;
}

inline void write_mirror_R_of_beta_csv(const std::filesystem::path& out, const mirrors_t& m) {
    std::ofstream f(out);
    f << "beta,R" << std::endl;
    for (const auto &entry : m.R_of_beta)
        f << std::to_string(entry.first) << ", " << std::to_string(entry.second) << std::endl;
}

inline void write_mirror_data_csv(const std::filesystem::path& out, const mirror_t& m, f_t step=.01) {
    std::ofstream f(out);
    f << "x1,y1,z1,x2,y2,z2" << std::endl;
    for (auto* p=&*m.begin();p+1<&*m.end();) {
        auto s=p+1;
        while (s+1<&*m.end() && glm::length(s->b-p->a)<step) ++s;
        f << serialize<8>(p->a) << "," << serialize<8>(s->b) << std::endl;
        p=s;
    }
}

}

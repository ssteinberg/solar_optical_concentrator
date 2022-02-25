
/*
 *  Copyright Shlomi Steinberg
 */

#pragma once

#include <atomic>
#include <mutex>

#include <glm/glm.hpp>

#include "aabb.hpp"
#include "ray.hpp"
#include "intersect.hpp"


namespace octree {

// P needs to conform to provide AABB() function
template <typename T, typename P, int max_level=10>
class octree {
public:
    using aabb_t = aabb<T>;

private:
    struct node_t {
        node_t *nnn{},*nnp{},*npn{},*pnn{},*npp{},*pnp{},*ppn{},*ppp{};
        std::vector<const P*> objs;
        std::mutex m;

        ~node_t() noexcept {
            if (nnn) delete nnn;
            if (nnp) delete nnp;
            if (npn) delete npn;
            if (pnn) delete pnn;
            if (ppn) delete ppn;
            if (pnp) delete pnp;
            if (npp) delete npp;
            if (ppp) delete ppp;
        }
    };

private:
    node_t root;
    aabb_t world;
    
protected:
    static std::size_t test(const aabb_t& aabb, const aabb_t& nbb, 
                            aabb_t &nnnbb, aabb_t &nnpbb, aabb_t &npnbb, aabb_t &pnnbb, aabb_t &nppbb, aabb_t &pnpbb, aabb_t &ppnbb, aabb_t &pppbb,
                            bool  &innn,   bool  &innp,   bool  &inpn,   bool  &ipnn,   bool  &inpp,   bool  &ipnp,   bool  &ippn,   bool  &ippp) noexcept {
        nnnbb = splice(nbb,{false,false,false}); nnpbb = splice(nbb,{false,false,true }); npnbb = splice(nbb,{false,true ,false}); pnnbb = splice(nbb,{true ,false,false});
        nppbb = splice(nbb,{false,true ,true }); pnpbb = splice(nbb,{true ,false,true }); ppnbb = splice(nbb,{true ,true ,false}); pppbb = splice(nbb,{true ,true ,true });
        innn = intersect(nnnbb,aabb); innp = intersect(nnpbb,aabb); inpn = intersect(npnbb,aabb); ipnn = intersect(pnnbb,aabb); 
        inpp = intersect(nppbb,aabb); ipnp = intersect(pnpbb,aabb); ippn = intersect(ppnbb,aabb); ippp = intersect(pppbb,aabb);
        return (!!innn)+(!!innp)+(!!inpn)+(!!ipnn)+(!!inpp)+(!!ipnp)+(!!ippn)+(!!ippp);
    }

    void create_node_and_insert(const P* obj, const aabb_t& aabb, const aabb_t& nbb, node_t* &n,std::size_t level) noexcept {
        std::atomic_thread_fence(std::memory_order_acquire);
        auto *expected = n;
        if (!expected) {
            auto *nn = new node_t;
            if (std::atomic_compare_exchange_strong_explicit(reinterpret_cast<std::atomic<node_t*>*>(&n),
                                                             &expected,
                                                             nn,
                                                             std::memory_order_release,
                                                             std::memory_order_relaxed))
                expected = nn;
            else
                delete nn;
        }
        return insert(obj,aabb,nbb,expected,level);
    }
    void insert(const P* obj, const aabb_t& aabb, const aabb_t& nbb, node_t *n,std::size_t level) noexcept {
        assert(aabb<=nbb);
        
        aabb_t nnnbb, nnpbb, npnbb, pnnbb, nppbb, pnpbb, ppnbb, pppbb;
        bool innn, innp, inpn, ipnn, inpp, ipnp, ippn, ippp;
        const auto its = test(aabb, nbb, nnnbb, nnpbb, npnbb, pnnbb, nppbb, pnpbb, ppnbb, pppbb, innn, innp, inpn, ipnn, inpp, ipnp, ippn, ippp);

        assert(its>0);
        if (its>1 || level==max_level) {
            std::unique_lock<std::mutex> l(n->m);
            n->objs.emplace_back(obj);
            return;
        }

             if (innn) create_node_and_insert(obj,aabb,nnnbb,n->nnn,level+1);
        else if (innp) create_node_and_insert(obj,aabb,nnpbb,n->nnp,level+1);
        else if (inpn) create_node_and_insert(obj,aabb,npnbb,n->npn,level+1);
        else if (ipnn) create_node_and_insert(obj,aabb,pnnbb,n->pnn,level+1);
        else if (inpp) create_node_and_insert(obj,aabb,nppbb,n->npp,level+1);
        else if (ipnp) create_node_and_insert(obj,aabb,pnpbb,n->pnp,level+1);
        else if (ippn) create_node_and_insert(obj,aabb,ppnbb,n->ppn,level+1);
        else if (ippp) create_node_and_insert(obj,aabb,pppbb,n->ppp,level+1);
    }
    bool erase(const P* obj, const aabb_t& aabb, const aabb_t& nbb, node_t *n,std::size_t level) noexcept {
        assert(aabb<=nbb);
        
        aabb_t nnnbb, nnpbb, npnbb, pnnbb, nppbb, pnpbb, ppnbb, pppbb;
        bool innn, innp, inpn, ipnn, inpp, ipnp, ippn, ippp;
        const auto its = test(aabb, nbb, nnnbb, nnpbb, npnbb, pnnbb, nppbb, pnpbb, ppnbb, pppbb, innn, innp, inpn, ipnn, inpp, ipnp, ippn, ippp);

        assert(its>0);
        if (its>1 || level==max_level) {
            std::unique_lock<std::mutex> l(n->m);
            for (auto it=n->objs.begin();it!=n->objs.end();++it) {
                if (*it == obj) {
                    n->objs.erase(it);
                    return true;
                }
            }
            return false;
        }

        bool removed = false;
             if (n->nnn && innn) removed = erase(obj,aabb,nnnbb,n->nnn,level+1);
        else if (n->nnp && innp) removed = erase(obj,aabb,nnpbb,n->nnp,level+1);
        else if (n->npn && inpn) removed = erase(obj,aabb,npnbb,n->npn,level+1);
        else if (n->pnn && ipnn) removed = erase(obj,aabb,pnnbb,n->pnn,level+1);
        else if (n->npp && inpp) removed = erase(obj,aabb,nppbb,n->npp,level+1);
        else if (n->pnp && ipnp) removed = erase(obj,aabb,pnpbb,n->pnp,level+1);
        else if (n->ppn && ippn) removed = erase(obj,aabb,ppnbb,n->ppn,level+1);
        else if (n->ppp && ippp) removed = erase(obj,aabb,pppbb,n->ppp,level+1);
        return removed;
    }

    template <typename Fun>
    void trace(const ray_t<T> &r, const Fun &f, const aabb_t& nbb, const node_t *n) const noexcept {
        for (const auto* obj : n->objs) {
            if (!f(obj))
                return;
        }

        aabb_t nnnbb, nnpbb, npnbb, pnnbb, nppbb, pnpbb, ppnbb, pppbb;
        bool innn, innp, inpn, ipnn, inpp, ipnp, ippn, ippp;
        test(aabb_t{}, nbb, nnnbb, nnpbb, npnbb, pnnbb, nppbb, pnpbb, ppnbb, pppbb, innn, innp, inpn, ipnn, inpp, ipnp, ippn, ippp);

        if (n->nnn && intersect(r,nnnbb)) trace(r,f,nnnbb,n->nnn);
        if (n->nnp && intersect(r,nnpbb)) trace(r,f,nnpbb,n->nnp);
        if (n->npn && intersect(r,npnbb)) trace(r,f,npnbb,n->npn);
        if (n->pnn && intersect(r,pnnbb)) trace(r,f,pnnbb,n->pnn);
        if (n->npp && intersect(r,nppbb)) trace(r,f,nppbb,n->npp);
        if (n->pnp && intersect(r,pnpbb)) trace(r,f,pnpbb,n->pnp);
        if (n->ppn && intersect(r,ppnbb)) trace(r,f,ppnbb,n->ppn);
        if (n->ppp && intersect(r,pppbb)) trace(r,f,pppbb,n->ppp);
    }

public:
    octree(const aabb_t& world) noexcept : world(world) {}

    void insert(const P* obj) noexcept {
        const auto &aabb = obj->AABB();
        assert(glm::all(glm::lessThanEqual(glm::vec<3,T>{0,0,0}, aabb.b-aabb.a)));
        insert(obj, aabb, world, &root, 0);
    }
    void erase(const P* obj) noexcept {
        if (!erase(obj, obj->AABB(), world, &root, 0))
            assert(false && "erase() failed: not found");
    }

    template <typename Fun>
    void trace(const ray_t<T> &r, const Fun &f) const noexcept {
        if (intersect(r,world))
            trace(r,f,world,&root);
    }
};

}

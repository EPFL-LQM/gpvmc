#ifndef _STATE_1_H
#define _STATE_1_H

#include <vector>
#include <iostream>
#include "defs.h"

class State_1 {
    protected:
        /*! \brief State in the enumerated single particle state basis.
         * m_part[f][p] is the index of the state p'th particle of flavour f
         * is in. The flavour is any good quantum number of the underlying mean-field
         * Hamiltonian such that is restricts the matrix of single particle overlap
         * one should consider when building the trial state. Written in this basis,
         * the state must be antisymmetrized in order to be really a state of fermions.
         */
        std::vector<uint_vec_t> m_part;
        /*! \brief State written in Fock space.
         * m_fock[f][s] is the index of the enumerated particle of flavour f in state s,
         * or m_part[f].size() if no flavour f particle occupies the s state. This is
         * redundant with m_part but is useful to have as well.
         */
        std::vector<uint_vec_t> m_fock;
        /*! \brief Size of the Fock space for each flavour f.
         * m_Nfs[f] is the number of different states a flavour f particle can be in.
         * It is given by m_Nfs[f]=m_fock[f].size()
         */
        uint_vec_t m_Nfs;
        /*! \brief Number of particles in flavour f.
         * It is given by m_Npt[f]=m_part[f].size()
         */
        uint_vec_t m_Npt;
        /*! \brief Fermionic sign
         * Sign of the state with respect to the initial one. The sign is updated upon
         * performing hopping.
         */
        int m_sign;

        /*! \brief Number of conserved flavours.
         * It is given by, for instance, m_Nfl=m_Nfs.size()
         */
        std::size_t m_Nfl;

        /*! \brief Utility function used in constructor of this class and derived class.
         */
        void build_base(const uint_vec_t& Npt,
                        const uint_vec_t& Nfs);

        /*! \brief Defines a sequence of hoppings to go from Fock state src to Fock state dest.
         */
        void get_hop_path(const uint_vec_t& src,
                          const uint_vec_t& dest,
                          const std::size_t& Npt,
                          hop_path_t& hop) const;
        /*! \brief Execute the hoppings defined in hop on Fock state src and modifies dest.
         */
        void get_hop_state(const hop_path_t& hop,
                           const uint_vec_t& src,
                           const std::size_t& Npt,
                           uint_vec_t& dest) const;

    public:
        /*! \brief Default constructor.
         */
        State_1();

        virtual ~State_1();

        /*! \brief Initialize the state with the particle state st.
         * st[f][p] is the state the f'th flavour
         * p'th particle is in.
         */
        void InitPart(const std::vector<uint_vec_t>& st);
        /*! \brief Initialize the state with the fock state st.
         * st[f][s]<m_Npt[f] if state s of flavour f is occupied.
         */
        void InitFock(const std::vector<uint_vec_t>& st);

        const uint_vec_t& GetNfs() const;
        const uint_vec_t& GetNpt() const;
        size_t GetNfl() const;
        const std::vector<uint_vec_t>& Getfs() const;
        const std::vector<uint_vec_t>& Getpt() const;
        int GetSign() const;

        /*! \brief Perform a sequence of hoppings
         */
        void Hop(const std::vector<hop_path_t>& hops);

        /*! \brief Calculate the fermionic sign change of a sequence of hoppings.
         */
        int HopSign(const std::vector<hop_path_t>& hops) const;

        /*! \brief Easy output
         */
        friend std::ostream& operator<<(std::ostream& left,
                                        const State_1& right);
};

#endif//_STATE_1_H

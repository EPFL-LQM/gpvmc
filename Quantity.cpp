#ifdef USEMPI
#include <mpi.h>
#endif
#include "Quantity.h"

Quantity::Quantity(const Stepper* stepper, FileManager* fm, std::string basename)
    :m_comm_size(1),m_comm_rank(0),m_meas(0),m_stepper(stepper),m_fm(fm),
     m_basename(basename)
{
#ifdef USEMPI
    MPI_Comm_size(MPI_COMM_WORLD,&m_comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&m_comm_rank);
#endif
}

Quantity::~Quantity()
{}

void Quantity::init()
{
    m_meas=0;
}

void Quantity::measure()
{
    m_meas++;
}

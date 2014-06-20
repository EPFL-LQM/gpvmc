#include "JastrowTrack.h"
#include "Stepper.h"
#include "FileManager.h"
#include "Jastrow.h"

JastrowTrack::JastrowTrack(const Stepper* stepper, FileManager* fm)
    :Quantity(stepper,fm,"JastrowTrack")
{}
JastrowTrack::~JastrowTrack() {}
void JastrowTrack::init()
{
    Quantity::init();
    m_vals.clear();
}
void JastrowTrack::save() const
{
    vector<double> amps(m_vals.size());
    BigDouble mean=0;
    for(size_t p=0;p<m_vals.size();++p)
        mean+=m_vals[p];
    mean/=BigDouble(double(m_vals.size()));
    for(size_t p=0;p<m_vals.size();++p)
        amps[p]=double(m_vals[p]/mean);
    m_fm->FileStream(m_basename).Write(m_vals.size(),1,amps.data());
}
void JastrowTrack::measure()
{
    m_vals.push_back(m_stepper->GetJas()->Jas());
}

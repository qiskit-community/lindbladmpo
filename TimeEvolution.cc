#include "TimeEvolution.h"
//____________________________________________________________________
void TimeEvolver::init(double tau, const AutoMPO &auto_L, Args args, int ord)
{
    order = ord;
    if (order > 4 || order < 2)
        cerr << "Error, TrotterOrder=" << order << " not implemented.\n", exit(1);

    argsApplyMPOtoRho = args; //Take the options given
    argsApplyMPOtoRho.add("Method", "Fit");
    //argsApplyMPOtoRho.add("Method", "DensityMatrix"); //Alternative method/algorithm to apply an MPO to an MPS. More precise.
    argsApplyMPOtoRho.add("Normalize", false);

    Cplx t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
    if (order == 2)
    {
        //Approx. with error O(tau^3)
        t1 = 0.5 * (1 + Cplx_i) * tau;
        t2 = 0.5 * (-1 + Cplx_i) * tau;
        t3 = 0;
        t4 = 0;
        expL1 = toExpH(auto_L, t1);
        expL2 = toExpH(auto_L, t2);
    }
    if (order == 3)
    {
        //Approx. with error O(tau^4). See Phys. Rev. B 96, 195117 (2017) appendix
        //Re(t1)= 1/4-1/12*3^(1/2)
        //Im(t1)=-1/4-1/12*3^(1/2)
        t1 = .10566243270259355887 - .39433756729740644113 * Cplx_i;
        t2 = Cplx_i * t1;
        t3 = conj(t2);
        t4 = Cplx_i * t3;
        t1 *= Cplx_i * tau;
        t2 *= Cplx_i * tau;
        t3 *= Cplx_i * tau;
        t4 *= Cplx_i * tau;
        expL1 = toExpH(auto_L, t1);
        expL2 = toExpH(auto_L, t2);
        expL3 = toExpH(auto_L, t3);
        expL4 = toExpH(auto_L, t4);
    }
    if (order == 4)
    {
        //Approx. with error O(tau^5). See Phys. Rev. B 96, 195117 (2017) appendix
        t1 = 0.25885339861091821723 + 0.04475613401114190287 * Cplx_i;
        t2 = -0.03154685814880379274 + 0.24911905427556321757 * Cplx_i;
        t3 = 0.19082905211066719664 - 0.23185374923210605447 * Cplx_i;
        t4 = 0.1637288148544367438753;
        t5 = conj(t3);
        t6 = conj(t2);
        t7 = conj(t1);
        t1 *= Cplx_i * tau;
        t2 *= Cplx_i * tau;
        t3 *= Cplx_i * tau;
        t4 *= Cplx_i * tau;
        t5 *= Cplx_i * tau;
        t6 *= Cplx_i * tau;
        t7 *= Cplx_i * tau;
        expL1 = toExpH(auto_L, t1);
        expL2 = toExpH(auto_L, t2);
        expL3 = toExpH(auto_L, t3);
        expL4 = toExpH(auto_L, t4);
        expL5 = toExpH(auto_L, t5);
        expL6 = toExpH(auto_L, t6);
        expL7 = toExpH(auto_L, t7);
    }
}
//____________________________________________________________________
void TimeEvolver::evolve(MPS &rho) const
{
    rho = applyMPO(expL1, rho, argsApplyMPOtoRho);
    rho.noPrime("Site");
    rho = applyMPO(expL2, rho, argsApplyMPOtoRho);
    rho.noPrime("Site");
    if (order == 3)
    {
        rho = applyMPO(expL3, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
        rho = applyMPO(expL4, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
    }
    if (order == 4)
    {
        rho = applyMPO(expL3, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
        rho = applyMPO(expL4, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
        rho = applyMPO(expL5, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
        rho = applyMPO(expL6, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
        rho = applyMPO(expL7, rho, argsApplyMPOtoRho);
        rho.noPrime("Site");
    }
}
//____________________________________________________________________
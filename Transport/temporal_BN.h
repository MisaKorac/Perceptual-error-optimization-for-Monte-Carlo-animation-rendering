//
//

#ifndef SLICEDOPTIM_NCUBESLICEDOPTIMALTRANSPORT_H
#define SLICEDOPTIM_NCUBESLICEDOPTIMALTRANSPORT_H

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <cstring>
#include "../Math/VecX.h"
#include "../Math/myMath.h"
#include "../Tools/iopointset.h"
#include "../Tools/my_utility.h"

int tileSize = 64;
int nb_frames = 1;
double sigma = 2.1;//1.9;//2.1
double sigma_temporal = 2.1;//1.9;//2.1
int kernel_size = 5;
int kernel_size_temporal = 5;
double mul = 3;
double r = 2;

double Rafal_weights[15] = {0.0395582850655909,0.398409765599820,0.404689211439515,0.177879764807883,0.0419136195197315,-0.00678276096076438,-0.0170052838869954,-0.0150087151111371,-0.0106941128889411,-0.00697594351554300,-0.00435253914170693,-0.00264862667445315,-0.00158680181327903,-0.000939965263123874,0.00354410282340402};

//double weight_temporal(int dt){
//    float epsilon = 0.0001;
//    //eturn std::exp(std::pow(std::log(dt*0.016+epsilon)-std::log(0.06),2)/(2*std::pow(sigma_temporal,2.0)));
//    return std::exp(-(dt*dt)/std::pow(sigma_temporal,2.0));
//}
double weight_temporal(int dt){
    float alpha = 0.2;//0.1
    //return std::pow(1-alpha,dt);
    return Rafal_weights[dt];
    //return std::pow(1-alpha,dt)*Rafal_weights[dt];
    //return std::exp(-(dt*dt)/std::pow(sigma_temporal,2.0));
    //return std::exp(-(dt*dt)/std::pow(sigma_temporal,2.0))*std::pow(1-alpha,dt);
}

double weight(int di,int dj, int dt){
    return std::exp(-(di*di+dj*dj)/std::pow(sigma,2.0))* weight_temporal(dt);
}

template <class VECTYPE>
void project(const std::vector<VECTYPE>& points,
            std::vector<std::pair<double, int>>& pointsProject,
            const VECTYPE& dir, VECTYPE& center, int functionSelection,
            double w_ref, int dim){

    int spp = points.size()/(tileSize*tileSize*nb_frames);
    // Weighting Function Selection
    int k = functionSelection%(tileSize*tileSize);
    int u = k/int(tileSize);
    int v = k%int(tileSize);
    int t = functionSelection/(tileSize*tileSize);
    
    if(tileSize == 1){
        kernel_size = 1;
    }

    float W_temporal =  ((float)rand() / RAND_MAX);

    for(int dt = -kernel_size_temporal+1; dt < 1; ++dt){
        for(int di = -kernel_size+1; di < kernel_size; ++di){
            for(int dj = -kernel_size+1; dj < kernel_size; ++dj){
                int x = ((u+di + tileSize)%tileSize);
                int y = ((v+dj + tileSize)%tileSize);
                int z = ((t+dt + nb_frames)%nb_frames);
                
                int pixel_indice = ((z*tileSize*tileSize + x*tileSize + y)*spp);

                if (weight(di,dj,std::abs(dt)) >= w_ref) {//  and weight_temporal(dt) >= W_temporal
                    for(int s = 0; s < spp; ++s){
                        VECTYPE p = points[pixel_indice+s];
                        double proj = toroidal_minus(p,center) * dir;

                        int indices = pixel_indice+s;

                        std::pair<double, int> val_indice = std::pair<double, int>(proj, indices);
                        pointsProject.push_back(val_indice);
                    }
                }
            }
        }
    }
}

/**
 *  \brief Get Optimal position at which to place \p nbSamples 1D sample to minimize OT cost to the Radon transform of a \p D dimensional ball
 *
 *  A common result in optimal transport is that 1D placement comes from the inverse of the CDF of target distribution
 * @param nbSamples Number of samples to compute position for
 * @param pos Output buffer in which to put the positions
 */
template <class VECTYPE>
inline void getInverseRadonNCube(int D, int nbSamples, std::vector<double>& pos,VECTYPE dir, VECTYPE& center){

    pos.resize(nbSamples);
    std::vector<double> point_projected;
    for (int i = 0; i < nbSamples * mul; i++)
    {
        VECTYPE p;
        for (int d = 0; d < D; d++)
        {
            p[d] = ((float)rand() / RAND_MAX);
        }
        
        point_projected.push_back(toroidal_minus(p ,center)* dir);        
    }

    std::sort(point_projected.begin(),point_projected.end());
    for (int i = 0; i < nbSamples; i++)
    {
        pos[i] = point_projected[i*mul + int(mul/2)];
    }
}

/**
 * Choose \p m function of weighting.
 *
 * @param function_slices Table of function selection to output.
 * @param m Number of directions to pick
 */
inline void chooseFunctionSlices(std::vector<int>& function_slices, int m){
    int offset = rand()%(tileSize*tileSize);
    for (int k = 0; k < m; ++k){
        int a = (double(k)/(m)) * (tileSize*tileSize) + offset;
        function_slices[k] = a%(tileSize*tileSize);
        function_slices[k] = rand()%(tileSize*tileSize*nb_frames);
    }
}

/**
 * Choose \p m function of weighting.
 *
 * @param function_slices Table of function selection to output.
 * @param m Number of directions to pick
 */
template <class VECTYPE>
inline void chooseWeightSlices(std::vector<double>& weight_slices, int m,VECTYPE c){
    for (int k = 0; k < m; ++k){
        double a = ((float)rand() / RAND_MAX)*weight(0,0,0);
        weight_slices[k] = a;
    }
}

/**
 * Compute optimal transport in 1D for direction \f$ \theta \f$ and \f$d_{j}\f$ being the 1D displacement of \f$\x^j\f$
 * that minimize the 1D sliced optimal transport along \f$ \theta \f$.
 *
 * Denoting $\sigma$ the permutations of the indices \f$\{j\}_{j=1..N}\f$ such that
 * \f$\bigl(\x^{\sigma(j)} \cdot \theta \bigr)_j\f$, is a sorted sequence of increasing values,
 * one can compute \f$d_{j}\f$ via
 * \f$ d_{j} = C_{\theta}^{-1}\left(\frac{\sigma(j)-\frac12}{N}\right)\,. \vspace*{-1mm}\f$
 *
 * @param dir Direction \f$ \theta \f$
 * @param points Table containing the points \f$ x_j \f$
 * @param shift Output the 1D shift to apply to minimize transport cost
 */
template<class VECTYPE>
inline void slicedStepNCube(const VECTYPE& dir, int functionSelection,
                            double W_ref, const std::vector<VECTYPE>& points,
                            std::vector<double>& shift)
{
    
    int N = points.front().dim();
    VECTYPE center;
    for (int i = 0; i < N; ++i){
        center[i] = 1.0*((float)rand() / RAND_MAX);
    }
    for (size_t i = 0; i < points.size(); ++i){
        shift[i] = 0.0f;
    }
    std::vector<std::pair<double, int>> pointsProject;
    project(points, pointsProject, dir, center, functionSelection, W_ref, N);
    if(pointsProject.size() >= 1){
        //Compute optimal 1D position for given number of points;
        std::vector<double> pos(pointsProject.size());
        getInverseRadonNCube(N, pointsProject.size(), pos, dir, center);
        std::sort(pointsProject.begin(), pointsProject.end(), [](const std::pair<double, int> &x,
                                        const std::pair<double, int> &y)
        {
            return x.first < y.first;
        });
        double grad_normalization = float(pointsProject.size())/points.size();
        double normalise_factor = 1.0;
        //Computes required shift to optimize 1D optimal transport
        for (size_t i = 1; i < pointsProject.size()-1; ++i) {
            //Compute shifting
            double s = pos[i] - pointsProject[i].first;
            normalise_factor = ((pos[i + 1] - pos[i - 1]) / 2.0) * pointsProject.size();
            shift[pointsProject[i].second] += grad_normalization*(s/normalise_factor);
        }
        double s = pos[0] - pointsProject[0].first;
        normalise_factor = (pos[1] - pos[0]) * pointsProject.size();
        shift[pointsProject[0].second] += grad_normalization*(s / normalise_factor);

        s = pos[pointsProject.size() - 1] - pointsProject[pointsProject.size() - 1].first;
        normalise_factor = (pos[pointsProject.size() - 1] - pos[pointsProject.size() - 1 - 1]) * pointsProject.size();
        shift[pointsProject[pointsProject.size() - 1].second] += grad_normalization*(s / normalise_factor);
    }
}

/**
 * Compute optimal transport in 1D for the \p directions and displace \f$ x_j \f$ by
 * \f$\pmb{\delta}^j \EqDef \frac{1}{K} \sum_{i=1}^K d_{i,j}\, \theta_i \vspace*{-1mm}\f$ with
 * with \f$ d_{i,j} \f$ being the displacement of \f$\x^j\f$ that minimize the 1D sliced optimal
 * transport along direction \f$ \theta_i \f$
 *
 * @param pointsOut Table containing the points \f$ x_j \f$
 * @param directions Table containing the \f$\theta_i\f$
 * @param shift Used to avoid having to allocate uge chunks of memory. Must be a vector of size m containing vectors of same size as \p pointsOut.
 * @param finalShift Used to avoid having to allocate huge chunks of memory Must be a vector of same size as \p pointsOut.
 * @return the Wasserstein cost of the current iteration.
 */
template <class VECTYPE>
inline void slicedOptimalTransportBatchCube(std::vector<VECTYPE>& pointsOut,
                                 const std::vector<VECTYPE>& directions,
                                 const std::vector<int>& function_slices,
                                 const std::vector<double>& weight_slices,
                                 std::vector<std::vector<double>>& shift,
                                 std::vector<VECTYPE>& finalShift,
                                 std::vector<VECTYPE>& m_adam,
                                 std::vector<VECTYPE>& v_adam,
                                 int t)
 {
    int m = directions.size();
    int nbPoints = pointsOut.size();
    int dim = pointsOut.front().dim();
    //Compute the shift along each direction
#pragma omp parallel for shared(directions, shift)
    for (int k = 0; k < m; ++k){
        for(double& v : shift[k]){
			v = 0.;
        }
        const VECTYPE& dir = directions[k];

        slicedStepNCube(dir,function_slices[k], weight_slices[k], pointsOut, shift[k]);
    }
        //Accumulate shift from all directions
#pragma omp parallel for
    for (int i = 0; i < nbPoints; ++i) {
        VECTYPE sh(finalShift[i].dim());
        memset(&sh[0], 0, finalShift[i].dim() * sizeof(sh[0]));
        for (int k = 0; k < m; ++k) {
            sh += shift[k][i] * directions[k];
        }

        for (int k = 0; k < dim; ++k) {
            double grad = (sh[k]/m)*-1;
            double lr = 0.5;//12.6;//8.6, 0.005
            double B1 = 0.9;
            double B2 = 0.99;
            double epsilon = 0.0001;
            m_adam[i][k] = B1*m_adam[i][k]+(1-B1)*grad;
            v_adam[i][k] = B2*v_adam[i][k]+(1-B2)*std::pow(grad,2);
            double m_hat = m_adam[i][k]/(1-std::pow(B1,t));
            double v_hat = v_adam[i][k]/(1-std::pow(B2,t));
            finalShift[i][k] = lr*m_hat/(std::sqrt(v_hat)+epsilon);
        }
    }
    //Displace points according to accumulated shift
#pragma omp parallel for
    for (int i = 0; i < nbPoints; ++i) {
        pointsOut[i] -= finalShift[i];

        for(int d = 0; d<finalShift[i].dim(); ++d){
            while(pointsOut[i][d]<0){
                pointsOut[i][d]+=1.0;
            }
            while(pointsOut[i][d]>1){
                pointsOut[i][d]-=1.0;
            }
        }
    }

}

void print_progress(double ratio){
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * ratio;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(ratio * 100.0) << " %\r";
    std::cout.flush();
}

/**
 *  \brief Computes an optimized point set to uniformly sample the unit N-Ball using sliced optimal transport
 *
 * @param pointsIn Contains input ND points
 * @param pointsOut Contains optimized ND points
 * @param nbIter Number of iterations
 * @param m Number of slice per iteration
 * @param seed random seed
 * @return the Sliced Wasserstein distance between the samples and the uniform distribution.
 */

template <class VECTYPE>
inline void slicedOptimalTransportNCube(const std::vector<VECTYPE>& pointsIn,
                                        std::vector<VECTYPE>& pointsOut,
                                        int nbIter,
                                        int m,
                                        int seed,
                                        int tileS,
                                        bool silent,
                                        int frames)
{
    tileSize = tileS;
    std::cout << frames << std::endl;
    nb_frames = frames;
    int N = pointsIn.front().dim();
    pointsOut = pointsIn;

    //Accumulation shift to be applied later on
    std::vector<std::vector<double>> shift(m, std::vector<double>(pointsOut.size()));
    std::vector<VECTYPE> m_adam(pointsOut.size(), VECTYPE(N));
    std::vector<VECTYPE> v_adam(pointsOut.size(), VECTYPE(N));
    std::vector<VECTYPE> finalShift(pointsOut.size(), VECTYPE(N));

    std::vector<VECTYPE> directions(m, VECTYPE(N));
    std::vector<int> function_slices(m);
    std::vector<double> weight_slices(m);
    
    srand (seed);
    //Iterate 1D Optimal Transport
    for (int i = 0; i < nbIter; i += 1){
        if(!silent){
            print_progress(double(i)/nbIter);
        }
        VECTYPE c;
        chooseDirectionsND(directions, m, seed);
        chooseFunctionSlices(function_slices, m);
        chooseWeightSlices(weight_slices, m,c);

        slicedOptimalTransportBatchCube(pointsOut, directions, function_slices, weight_slices, shift, finalShift,m_adam,v_adam,i+1);
    }
    print_progress(1.0);

}


#endif //SLICEDOPTIM_NCUBESLICEDOPTIMALTRANSPORT_H

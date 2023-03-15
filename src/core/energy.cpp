/**
 * @file energy.cpp
 * @author Bohan Cao (2110313@mail.nankai.edu.cn)
 * @brief 
 * @version 0.0.1
 * @date 2023-03-08
 * 
 * @copyright Copyright (c) 2023
 */
#include "core/energy.h"

#include "core/math.h"

namespace CoreEnergy {

double mean_curvature(const CoreMath::Array<CoreMath::Vector>& others) {
  CoreMath::Vector H_vec;
  int size = 0;

  for (int i = 0; i < others.size(); i++) {
    int left = (i==0) ? others.size()-1 : i-1;
    int right = (i==others.size()-1) ? 0 : i+1;

    // vectors for a node and its related bonds
    CoreMath::Vector left_up = others[left] - others[i];
    CoreMath::Vector right_up = others[right] - others[i];

    // tan(acos(x)) = sqrt(1/x/x-1)
    double cos_angle1 = left_up*others[left]/CoreMath::mod(left_up)/
        CoreMath::mod(others[left]);
    double cos_angle2 = right_up*others[right]/CoreMath::mod(right_up)/
        CoreMath::mod(others[right]);
    H_vec += (1/Kokkos::sqrt(1/cos_angle1/cos_angle1-1) + 
        1/Kokkos::sqrt(1/cos_angle2/cos_angle2-1)) * others[i];
    size += CoreMath::mod(CoreMath::cross(others[i], others[right])) / 2;
  }

  return CoreMath::mod(H_vec)*3/4/size;
}

double gaussian_curvature(const CoreMath::Array<CoreMath::Vector>& others) {
  // Size around this node, total angle
  double angles = 2*PI, size = 0;
  for (int i=0; i<others.size(); i++) {
    int right = (i==others.size()-1) ? 0 : i+1;

    angles -= Kokkos::acos(others[right]*others[i]/CoreMath::mod(others[right])/
        CoreMath::mod(others[i]));
    size += CoreMath::mod(CoreMath::cross(others[i], others[right])) / 2;
  }

  return angles*3/size;
}
    
double curvature_energy(const double* para, const CoreMath::Array<CoreMath::Vector>& others) {
  CoreMath::Vector H_vec;
  double size = 0, angles = 2*PI;

  for (int i = 0; i < others.size(); i++) {
    int left = (i==0) ? others.size()-1 : i-1;
    int right = (i==others.size()-1) ? 0 : i+1;

    // vectors for a node and its related bonds
    CoreMath::Vector left_up = others[left] - others[i];
    CoreMath::Vector right_up = others[right] - others[i];

    // tan(acos(x)) = sqrt(1/x/x-1)
    double cos_angle1 = left_up*others[left]/CoreMath::mod(left_up)/
        CoreMath::mod(others[left]);
    double cos_angle2 = right_up*others[right]/CoreMath::mod(right_up)/
        CoreMath::mod(others[right]);
    H_vec += (1/Kokkos::sqrt(1/cos_angle1/cos_angle1-1) + 
        1/Kokkos::sqrt(1/cos_angle2/cos_angle2-1)) * others[i];
    size += CoreMath::mod(CoreMath::cross(others[i], others[right])) / 2;
    angles -= Kokkos::acos(others[right]*others[i]/CoreMath::mod(others[right])/
        CoreMath::mod(others[i]));
  }

  return para[0]*((H_vec*H_vec)*9/8/size/size - angles*3/size);
}

CoreMath::Array<CoreMath::Vector> curvature_gradient(const double* para, 
    const CoreMath::Array<CoreMath::Vector>& others) {
  // angle, size*2
  double angles = 2*PI, size2 = 0;
  CoreMath::Vector H_vec;

  CoreMath::Array<CoreMath::Vector> result(others.size());

  // compute a*b, |a\times b|, (|a||b|)^2, (c-a)*c, |(c-a)\times c|, |c-a||c|, 
  // (b-a)*b, |(b-a)\times b|, |b-a||b|
  auto s = others.size();
  CoreMath::Array<double> a0(s), a1(s), a2(s), b0(s), b1(s), b2(s), c0(s), c1(s), c2(s);
  for (int i=0; i<others.size(); i++) {
    int il = (i==0) ? others.size()-1 : i-1;
    int ir = (i==others.size()-1) ? 0 : i+1;

    // vectors for a node and its related bonds
    CoreMath::Vector il_i = others[il] - others[i];
    CoreMath::Vector ir_i = others[ir] - others[i];

    a0[i] = others[i]*others[ir];
    a1[i] = CoreMath::mod(CoreMath::cross(others[i], others[ir]));
    a2[i] = Kokkos::sqrt((others[ir]*others[ir])*(others[i]*others[i]));
    b0[i] = others[il]*il_i;
    b1[i] = CoreMath::mod(CoreMath::cross(others[il], il_i));
    b2[i] = Kokkos::sqrt((others[il]*others[il])*(il_i*il_i));
    c0[i] = others[ir]*ir_i;
    c1[i] = CoreMath::mod(CoreMath::cross(others[ir], ir_i));
    c2[i] = Kokkos::sqrt((others[ir]*others[ir])*(ir_i*ir_i));

    H_vec += (b0[i]/b1[i] + c0[i]/c1[i]) * others[i];
    angles -= Kokkos::acos(a0[i]/a2[i]);
    size2 += a1[i];
  }

  // Gaussian curvature and mean curvature
  double G = 6*angles/size2, H = CoreMath::mod(H_vec)*3/2/size2;
  double G_4H2 = G - 4*H*H;
  H_vec = H_vec / CoreMath::mod(H_vec);

  // compute gradient
  for (int i=0; i<others.size(); i++) {
    int il = (i==0) ? others.size()-1 : i-1;
    int ir = (i==others.size()-1) ? 0 : i+1;

    // vectors for a node and its related bonds
    CoreMath::Vector il_i = others[il] - others[i];
    CoreMath::Vector ir_i = others[ir] - others[i];

    double k = para[0]/a1[i]/size2;
    double k1 = 6 + G_4H2*a0[i];
    double k2 = 6*a0[i]/a2[i]/a2[i] + G_4H2;
    result[i] += k*(k2*(others[ir]*others[ir])*others[i] - k1*others[ir]);
    result[ir] += k*(k2*(others[i]*others[i])*others[ir] - k1*others[i]);

    k = 6*para[0]*H/size2;
    k1 = b2[i]*b2[i]/b1[i]/b1[i]/b1[i]*(others[i]*H_vec);
    k2 = c2[i]*c2[i]/c1[i]/c1[i]/c1[i]*(others[i]*H_vec);
    result[i] += k*((b0[i]/b1[i] + c0[i]/c1[i]) * H_vec -
        k1*(others[il] - b0[0]*il_i/(il_i*il_i)) -
        k2*(others[ir] - c0[0]*ir_i/(ir_i*ir_i)));
    result[il] += k*k1*((1-b0[0]/(others[il]*others[il])) * others[il] +
        (1-b0[0]/(il_i*il_i))*il_i);
    result[ir] += k*k2*((1-c0[0]/(others[ir]*others[ir])) * others[ir] +
        (1-c0[0]/(ir_i*ir_i))*ir_i);
  }

  return result;
}

} // namespace CoreEnergy
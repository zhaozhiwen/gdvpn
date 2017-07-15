Code for J/Psi production in deuteron breakup

==== how to compile, tested on CentOS6.5_x86_64 ========

git clone https://github.com/zhaozhiwen/gdvpn.git

wget http://www.feynarts.de/cuba/Cuba-4.2.tar.gz

tar zxf Cuba-4.2.tar.gz

cd Cuba-4.2

./configure

make

cd ../

ln -s Cuba-4.2/libcuba.a

ln -s Cuba-4.2/cuba.h

cmake ./

make

== how to run and result ================================

run like "./driver.elf -E9 -v"

result will be in text file "test" with crossection and final particles mom and angle

see more options and output format in "driver.cc"

== info =================================================

I'm writing to send you some code I wrote during graduate school for computing the cross section of incoherent J/Psi photoproduction, in the hope that it will be useful in a Monte Carlo event generator.

The central part of the code is gdvpn.cc, which has a routine to compute the five-fold differential cross section with respect to p_V, Omega_V, and Omega_p. It is a function of t, p_n, theta_nl, phi_nl, and phi_v, where l is the three-momentum transferred from the incident photon into the deuteron. There is also a method init_gdvpn which can be used to set the initial kinematics of the photon and deuteron, and can be used to turn on or off specific impulse or rescattering contributions to the cross section.

There is also a method vnset in vmeson_params.cc which can be used to set properties of the meson-nucleon interaction, such as which meson (it can do J/Psi or phi(1020)), the total nucleon-meson cross section, the slope parameter B in dsig/dt=Ae^{Bt}. There is also an integer ID for picking a model for describing meson photoproduction from the nucleon; for J/Psi, it's 2 for two-gluon exchange and 3 for three-gluon. If you set B to be negative, then the vector meson dominance model is used to relate photoproduction and VN scattering.

The code is organized in a strange way (and I wouldn't do this with something I wrote today); each .cc file is meant to be compiled into an object separately, and then linked to a driver program. Static variables and methods are effectively used as private variables and methods. In retrospect, I should have used classes for private members/methods instead. If you would like, I could refactor the code this way, so you don't have to compile a bunch of objects separately, but this would take me a bit of time to do.

This code uses the Cuba integration library, cf. http://www.feynarts.de/cuba/ , to do numerical integration, so this is a dependency.

I've included a CMakeLists.txt file too, so you can use cmake to build the libraries and a driver program. The driver program was meant to create plaintext files with tables of cross sections at given kinematics in order to make plots. It takes a few optional inputs, which are described in the comments in driver.cc.

== log ===================================================

based on paper https://journals.aps.org/prc/abstract/10.1103/PhysRevC.88.044604

first checkin code from Adam J. Freese <afreese@anl.gov>

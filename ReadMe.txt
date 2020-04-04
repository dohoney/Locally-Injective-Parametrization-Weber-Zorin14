******************************************
01/02/2016

Disclaimer: This is an unofficial implementation of the SIGGRAPH 2014 paper
"Locally Injective Parametrization with Arbitrary Fixed Boundaries" authored by
Ofir Weber and Denis Zorin. This is not the implementation used by the authors
to produce the results that appear in the paper. The implementation was done by
Moshiko Schwortz and Alon Bright in the framework of an undergraduate student
project in the Engineering department in Bar Ilan University. This implementation
does not include the extensions to variable metric and extremal quasiconformal
maps that appears in sections 5.2 and 5.3 of the paper. These can dramatically
improve the quality of the results, though the current implementation possesses
the theoretical guarantees on local injectivity.

The use of this application is limited to academic use only!

The code is provided as-is and without any guarantees.

For any questions or comments please contact Alon Bright ( alonbright@gmail.com ).

------------------------------------------
A Visual Studio 2013 (MSVC 12.0) project is provided for easy compilation on
Windows machines though the code should be (not verified) platform
independent.

We used the following prerequisites: GMM, MATLAB R2013b, Boost 1_59,
CGAL 4.7

------------------------------------------
Environment Variables that needs to be set:

GMM_INCLUDE_DIR ----- (%your GMM folder path%)\gmm-4.2\include

MATLAB_64_DIR ------- C:\Program Files\MATLAB\R2013b

CGAL_DIR ------------ C:\Program Files\CGAL-4.7

BOOST_INCLUDEDIR ---- C:\boost\boost_1_59_0

BOOST_LIBRARYDIR ---- C:\boost\boost_1_59_0\lib64-msvc-12.0
------------------------------------------

In addition, you need to extract the folder "matlabScripts" and add it to your MATLAB path.

Tutorial.pdf contains further instructions on how to use the application.

 ******************************************

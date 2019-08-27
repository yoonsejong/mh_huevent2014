 A MATLAB implementation Psychology and Art Theory-based Features (MH)
 
 This set of MATLAB scripts extract the MH features [1] used in [2]. This code does NOT run itself out-of-the-box since it depends on many other programs and the AIC dataset. While we tried our best to make the job as simplified as possible, some of them must be downloaded and compiled separately due to license/compatibility issues. The code to compute Tamura texture features were ported from a C code of [7] under the MIT license included.

 You must cite [1] in order to use these features. In addition, please  consider citing [2] if you found this code useful.

 Please let me know if you found any bug(s).

 ==========
 HOW TO USE
 ==========

  STEP 1. Download hsy_rgb.zip from http://allan.hanbury.eu/.
          The file is located under Colour Resources (as of 2015/02/10).
          Copy the two files into the same directory of this file:
           hsy2rgb.m, rgb2hsy.m

  STEP 2. Download FastEMD solver from [3]. Compile. Copy following two
          files into the same directory of this file:
           emd_hat_gd_metric_mex.mex*, emd_hat_mex.mex* 

  STEP 3. Download Color Descriptor from [4]. Compile. Copy following 
          files into the same directory of this file:
           mexColorNaming.mex*, w2c.mat
 
  STEP 4. Download color space converter from [5]. Compile. Copy following
          files into the same directory of this file:
           colorcalc.mex*, colorspace.m

  STEP 5. Download RGB histogram calculator from [6]. Copy following file
          into the same directory of this file:
           rgbhist_fast.m

  STEP 6. You need the MATLAB Computer Vision Toolbox. If not, you need 
          to obtain a Viola-Johns Face Detector implementation somewhere.

  STEP 7. You need to obtain AIC dataset to properly train the Fuzzy
          C-Means. While I included pretrained models (AIC_fcm_*.mat), 
          it is strongly recommended to train them again. Particularly,
          if you were to reproduce [1], you must re-train it using IAPS
          dataset instead of the pretrained models. Please check out
          the codes fuzzyFuncLearn*.m for more details. 

  =======
  LICENSE
  =======
    A MATLAB implementation Psychology and Art Theory-based Features (MH)
    Copyright (C) 2014 Sejong Yoon (sjyoon@cs.rutgers.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 =========
 REFERENCE
 =========
  [1] J. Machajdik and A. Hanbury,
      Affective Image Classification using Features Inspired by Psychology and Art Theory,
      ACM Intl. Conf. on Multimedia (MM), 2010.
      URL: http://www.imageemotion.org/
  [2] S. Yoon and V. Pavlovic, 
      Sentiment Flow for Video Interestingness Prediction,
      ACM Intl. Conf. on Multimedia (MM) Workshop (HuEvent), 2014.
  [3] O. Pele and M. Werman,
      Fast and Robust Earth Mover's Distances,
      ICCV 2009.
      URL: http://www.ariel.ac.il/sites/ofirpele/FastEMD/
  [4] J. van de Weijer, C. Schmid, J. Verbeek, D. Larlus
      Learning Color Names for Real-World Applications,
      IEEE Trans. in Img. Proc. (TIP), vol 18 (7):1512-1524, 2009.
      URL: http://cat.uab.es/~joost
  [5] P. Getreuer
      MATLAB Central File Exchange - Colorspace Transformations
      URL: .../28790-colorspace-transformations
  [6] M. K. Reddy,
      MATLAB Central File Exchange - Color Histogram of an RGB Image.
      URL: .../43630-color-histogram-of-an-rgb-image
  [7] T. Minka and R. W. Picard,
      A Sample Implementation of Tamura Texture Feature in C.
      URL: http://vismod.media.mit.edu/pub/tpminka/features/

Last Updated: February 10, 2015

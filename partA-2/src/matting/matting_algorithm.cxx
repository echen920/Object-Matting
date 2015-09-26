// 
//  USE THIS FILE FOR YOUR MATTING-SPECIFIC CODE
//
//  DO NOT MODIFY THIS FILE ANYWHERE EXCEPT WHERE 
//  EXPLICITLY NOTED!!!!
//




#include "matting.h"


//
//  The Triangulation Matting algorithm
//
//

bool matting::compute(void)
{
	if ((alpha_computed_ == true) && (outdated_ == false))
		// the results have already been computed, so
		// we have nothing to do
		return true;

	if (((bool) comp_1_ == false) || 
		((bool) comp_2_ == false) || 
		((bool) back_1_ == false) || 
		((bool) back_2_ == false))
		// we do not have enough information yet to run the
		// algorithm, we have nothing to do
		return false;

	// Ok, we have enough information to proceed
	
	// first let us allocate space for the result images
	alpha_.set_size(ni_, nj_);
	object_.set_size(ni_, nj_);


	//////////////////////////////////////////////////
	// PLACE YOUR CODE BETWEEN THESE LINES          //
	//////////////////////////////////////////////////

	int i,j;
	vil_image_view<vil_rgb<vxl_byte> > temp;
	// set object and alpha to the size of the input
	object_.set_size(comp_1_.ni(), comp_1_.nj());
	alpha_.set_size(comp_1_.ni(), comp_1_.nj());
	for (unsigned j = 0; j < comp_1_.nj(); ++j) {
		for (unsigned i = 0; i < comp_1_.ni(); ++i) {
				/* create a 6X4 table of the coefficients from the triangulation
				matting equation with their values.  */
				vnl_matrix<double> Tri(6,4,0.0);
				Tri(0,0) = 1.0;
				Tri(0,3) = -back_1_(i,j).r;
				Tri(1,1) = 1.0;
				Tri(1,3) = -back_1_(i,j).g;
				Tri(2,2) = 1.0;
				Tri(2,3) = -back_1_(i,j).b;
				Tri(3,0) = 1.0;
				Tri(3,3) = -back_2_(i,j).r;
				Tri(4,1) = 1.0;
				Tri(4,3) = -back_2_(i,j).g;
				Tri(5,2) = 1.0;
				Tri(5,3) = -back_2_(i,j).b;
				// create a table for each the inverse and transpose of Tri
				vnl_matrix<double> Tri_inv;
				vnl_matrix<double> Tri_trans;
				Tri_trans = Tri.transpose();
				// Pseudo Inverse = [(A^(T)A)^(-1)]*A^(T)
				Tri_inv = (vnl_matrix_inverse<double>(Tri_trans*Tri)*Tri_trans);
				vnl_matrix<double> comp(6,1,0.0);
				// create a 6X1 table for R,G,B delta in both comp1 and comp2
				comp(0,0) = comp_1_(i,j).r - back_1_(i,j).r;
				comp(1,0) = comp_1_(i,j).g - back_1_(i,j).g;
				comp(2,0) = comp_1_(i,j).b - back_1_(i,j).b;
				comp(3,0) = comp_2_(i,j).r - back_2_(i,j).r;
				comp(4,0) = comp_2_(i,j).g - back_2_(i,j).g;
				comp(5,0) = comp_2_(i,j).b - back_2_(i,j).b;

				/* object pixel color = composite pixel color * inverse
				of coefficient matrix. */
				vnl_matrix<double> res = Tri_inv*comp;
				// restrict alpha value to be between 0 and 1
				if (res(3,0) > 1.0) {
					res(3,0) = 1.0;
				} else if (res(3,0) < 0.0) {
					res(3,0) = 0.0;
				}
				// restrict pixel values to be between 0 and 255
				for (int z=0; z < 3; ++z) {
					if (res(z,0) > 255) {
						res(z,0) = 255;
					} else if (res(z,0) < 0) {
						res(z,0) = 0;
					}
				}
				object_(i,j).r = res.get(0,0) * res.get(3,0);
				object_(i,j).g = res.get(1,0) * res.get(3,0);
				object_(i,j).b = res.get(2,0) * res.get(3,0);
				alpha_(i,j) = res.get(3,0)*255.0;
			}
		}
	/////////////////////////////////////////////////


	alpha_computed_ = true;
	outdated_ = false;

	return true;
}

bool matting::compute_composite(vil_image_view<vil_rgb<vxl_byte> > input_im,
				vil_image_view<vil_rgb<vxl_byte> > &output_im)
{
	//////////////////////////////////////////////////
	// PLACE YOUR CODE BETWEEN THESE LINES          //
	//////////////////////////////////////////////////
	int i,j;
	double a;
	// set output image to be the size of the input image
	output_im.set_size(input_im.ni(), input_im.nj());
	// compute the output with the equation C = C_o*alpha + (1-alpha)*C_k
	for (unsigned j = 0; j < input_im.nj(); ++j) {
		for (unsigned i = 0; i < input_im.ni(); ++i) {
			a = alpha_(i,j)/255.0;
			output_im(i,j) = object_(i,j)*a + input_im(i,j)*(1-a);
		}
	}
	new_back_ = input_im;
	new_comp_ = output_im;
	new_composite_computed_ = true;
	//////////////////////////////////////////////////

	return true;
}

//////////////////////////////////////////////////
// PLACE ANY ADDITIONAL CODE BETWEEN THESE LINES//
//////////////////////////////////////////////////

//////////////////////////////////////////////////

#include<iostream>
using namespace std;

double func(double num){
	return 2*num +1;
}

double bisectorMethod(double left, double right, double err, double depth){
	if (depth == 0){
		return (left + right) / 2.0;
	}

	double Y_left, Y_right, Y_mid;
	double mid = (left + right) / 2.0;

	Y_left = func(left);
	Y_mid = func(mid);
	Y_right = func(right);

	if (Y_mid == 0) {
		return mid;
	}
	else if ((Y_mid + err > 0) && (Y_mid - err < 0))
	{
		return mid;
	}
	else{
		if (Y_left * Y_mid >= 0 ){
			return bisectorMethod(mid, right, err, depth - 1);
		}
		else{
			return bisectorMethod(left, mid, err, depth - 1);
		}
	}
}

	

int main(){
	double result, left, right, err, max_iterations = 10;
	
	cout << "Left Guess: \n";
	cin  >> left;

	cout << "Right Guess: \n";
	cin  >> right;

	cout << "Error Margin: \n";
	cin  >> err;
	
	result = bisectorMethod(left, right, err, max_iterations);
	
	cout << "A zero is located at " << result << endl << "\n";

	double trash;
	cin >> trash;
	return(0);
}

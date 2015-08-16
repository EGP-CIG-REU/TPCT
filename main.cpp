#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

// fstream is the file io

// explicitly state which variables to be used from the namespace std such as ifstream

using std::ifstream;
using std::cout;
using std::cin;
using std::string;

int main(int argc, char *argv[])
{
  //cout << argc << "\n";
	
	if (argc < 5)
	  return 1;
	//cout << argv[2] << "\n";
	
	// Output the number of arguments 'argc' passed to 'checker.exe'

	cout << argc << "\n";

	// we are going to compare two files at a time

	cout << (argc - 1) / 2 << "\n";
	
	// we are going to compare 'fileCount' pairs of files

	int fileCount = (argc - 1) / 2;

	int pointCount, fpointCount;

	// pointers to arrays of pointers, either a two dimensional array or ..

	double ***positions, ***positionsNext;
	int **ids, **idsNext;
	int **remapPoints;
	int **remapPointsNext, **mapPointsNext;

	// Initialize the double pointer arrays

	positions = new double**[fileCount];
	positionsNext = new double**[fileCount];

	// Take the outside pointer and make remappoints point to an array of integers.
    // to compenste for the fact that the ids are not in the same order in each file.

	remapPoints = new int*[fileCount];

	mapPointsNext = remapPointsNext = new int*[fileCount];

	ids = new int*[fileCount];
	idsNext = new int*[fileCount];
	
	// each pointer in myFiles points to the ifstream for the file. containing the initial position for a
	// given run while myFilesNext points to the file containing the positions of the
	// points at the end of that computation

	ifstream *myFiles = new ifstream[fileCount];
	ifstream *myFilesNext = new ifstream[fileCount];
	
	{
	  int i = 0;

		cout << "---First File!---\n";

		string testStr, testStr2;

		// Initialize the first two files in the current pair of files

		myFiles[i].open(argv[(2*i)+1]);
		myFilesNext[i].open(argv[(2*i)+2]);

		// Look for

		//Initial testStr to empty string then continues until :w

		for (testStr = ""; testStr != "Piece";  )
		{
			myFiles[i].ignore(1024, '<'); // Skip until we find a left angle brackets.
			myFiles[i] >> testStr;
		}

		for (testStr = ""; testStr != "Piece";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}

		// Use \" to indicate we're  looking for a " before the 52.

		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFiles[i].ignore(1024, ' '); // Skip empty space
			myFiles[i] >> testStr;        // testStr now contains NumberOfPoints="52" including the quotes
		}

		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFilesNext[i].ignore(1024, ' ');
			myFilesNext[i] >> testStr;
		}

		// put in the digits 52 into fpointCount, stop at the double quote

		fpointCount = stoi(testStr.substr(16,16),nullptr);

		//cout << stoi(testStr.substr(16,16),nullptr) << "\n";

		// set up arrays Positions , etc with fpointCount doubles

		positions[i] = new double*[fpointCount];
		positionsNext[i] = new double*[fpointCount];
		remapPoints[i] = new int[fpointCount];
		mapPointsNext[i] = new int[fpointCount];
		ids[i] = new int[fpointCount];
		idsNext[i] = new int[fpointCount];
		
		// Skip ahead in the file until the string is "DataArray"

		for (testStr = ""; testStr != "DataArray";  )
		{
			myFiles[i].ignore(1024, '<');
			myFiles[i] >> testStr;
		}

		for (testStr = ""; testStr != "DataArray";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		
		for (int j = 0; j < fpointCount; j++)
		{
		  remapPoints[i][j] = j;  // the jth point's id in the ith file


		  myFiles[i].ignore(1024, '\n'); // go the line after the line that begins with "<DataArray"
		  myFilesNext[i].ignore(1024, '\n');  // sip hte newline on the line positions[i][j]


		  // initialize the array "positions[i][j]" to an array of three doubles.

          positions[i][j] = new double[3];
		  positionsNext[i][j] = new double[3];


		  myFiles[i] >> positions[i][j][0];  // x-position
		  myFiles[i] >> positions[i][j][1];  // y-position
		  myFiles[i] >> positions[i][j][2];  // z-position

		  myFilesNext[i] >> positionsNext[i][j][0];
		  myFilesNext[i] >> positionsNext[i][j][1];
		  myFilesNext[i] >> positionsNext[i][j][2];

		  // output each position to the screen - just as a double check

			cout << " Position [" << i << "][" << j << 
			  "]: " << positions[i][j][0] <<
			  ", "  << positions[i][j][1] << 
			  ", "  << positions[i][j][2] << "\n";

			cout << " Next Position [" << i << "][" << j << // Next Position is actually the final position.
			  "]: " << positionsNext[i][j][0] <<
			  ", "  << positionsNext[i][j][1] << 
			  ", "  << positionsNext[i][j][2] << "\n";
		}
		
		// Now find the line containing the string "id"

		for (testStr = ""; testStr != "Name=\"id\"";  )
		{		 
			// We do this by looking for lines that contain the string "DataArray"

			for (testStr = ""; testStr != "DataArray";  )
			{
				myFiles[i].ignore(1024, '<');
				myFiles[i] >> testStr;
			}

			// for each line containing (starting with "<DataArray" find the string "Name" in that line

			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFiles[i].ignore(1024, ' ');
				myFiles[i] >> testStr;
			}

			// Go back to the top of the outer loop and check to see if the testStr is "Name="id"


		}
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{		 
			for (testStr = ""; testStr != "DataArray";  )
			{
				myFilesNext[i].ignore(1024, '<');
				myFilesNext[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFilesNext[i].ignore(1024, ' ');
				myFilesNext[i] >> testStr;
			}
		}
		for (int j = 0; j < fpointCount; j++)
		{
			myFiles[i].ignore(1024, '\n');
			myFilesNext[i].ignore(1024, '\n');

			// read in the ids for the pair of files

			myFiles[i] >> ids[i][j];
			myFilesNext[i] >> idsNext[i][j];

			// print them out

			cout << " id [" << i << "][" << j << 
			  "]: " << ids[i][j] <<
			  ", "  << idsNext[i][j] << "\n";
		}

		// check the ids in the first file in the pair (initial positions) against
		// the second file in the pair (final positions)

		for (int j = 0; j < fpointCount; j++)
		{
		  if (ids[i][j] == idsNext[i][j])

			  // In the ith file pair the j-th particle in the first file corresponds to the
			  // j-th point in the second file.

			  mapPointsNext[i][j] = j;
			else
			{
			  for (int k = 0; k < fpointCount; k++)
				{

				  if (ids[i][j] == idsNext[i][k])

					  // Otherwise in the ith file pair the j-th particle in the first file corresponds to the
					  // k-th point in the second file.

					  mapPointsNext[i][j] = k;
				}
			}
		}
		
		cout << "\n";
	}
	
	// Now do the same thing for the remaining pair of files

	for (int i = 1; i < fileCount; i++)
	{

		cout << "----File [" << i << "]!----\n";

		string testStr;

		myFiles[i].open(argv[(2*i)+1]);
		myFilesNext[i].open(argv[(2*i)+2]);

		for (testStr = ""; testStr != "Piece";  )
		{
			myFiles[i].ignore(1024, '<');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr != "Piece";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFiles[i].ignore(1024, ' ');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr.substr(0,16) != "NumberOfPoints=\"";  )
		{
			myFilesNext[i].ignore(1024, ' ');
			myFilesNext[i] >> testStr;
		}
		
		// cout << "Test1\n";
		
		pointCount = stoi(testStr.substr(16,16),nullptr);
		//cout << stoi(testStr.substr(16,16),nullptr) << "\n";
		positions[i] = new double*[pointCount];
		positionsNext[i] = new double*[pointCount];
		remapPoints[i] = new int[fpointCount];
		mapPointsNext[i] = new int[pointCount];
		ids[i] = new int[pointCount];
		idsNext[i] = new int[pointCount];
		
		for (int j = pointCount; j < fpointCount; j++)
		{
		  remapPoints[i][j] = -1;
		}
		
		for (testStr = ""; testStr != "DataArray";  )
		{
			myFiles[i].ignore(1024, '<');
			myFiles[i] >> testStr;
		}
		for (testStr = ""; testStr != "DataArray";  )
		{
			myFilesNext[i].ignore(1024, '<');
			myFilesNext[i] >> testStr;
		}
		// cout << "test2\n";

		// Same as loop for the first pair of files until this point in the new loop

		for (int j = 0; j < pointCount; j++)
		{
			myFiles[i].ignore(1024, '\n');
			myFilesNext[i].ignore(1024, '\n');

			positions[i][j] = new double[3];
			positionsNext[i][j] = new double[3];

			myFiles[i] >> positions[i][j][0];
			myFiles[i] >> positions[i][j][1];
			myFiles[i] >> positions[i][j][2];

			myFilesNext[i] >> positionsNext[i][j][0];
			myFilesNext[i] >> positionsNext[i][j][1];
			myFilesNext[i] >> positionsNext[i][j][2];
			
			// Now check to see if the initial position of each particle in the file myFiles[i] is the same
			// as the initial position of the jth particle in the myFiles[0], the very first initial position
			// file.

			if (j < fpointCount)
			{	
				if (positions[i][j][0] == positions[0][j][0] &&
						positions[i][j][1] == positions[0][j][1] &&
						positions[i][j][2] == positions[0][j][2])
				{
					remapPoints[i][j] = j;
				}
				else
				{
					//remapPoints[i][j] = -1;

					// if not then find the integer k such that a positions[i][j][0] == positions[0][k][0], etc
					for (int k = 0; k < fpointCount; k++)
					{
						if (positions[i][j][0] == positions[0][k][0] &&
								positions[i][j][1] == positions[0][k][1] &&
								positions[i][j][2] == positions[0][k][2])
						{
							if (k < fpointCount)
								remapPoints[i][k] = j;
						}
					}
				}
			}
			//cout << "test2.1\n";
		}
		
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{		 
			for (testStr = ""; testStr != "DataArray";  )
			{
				myFiles[i].ignore(1024, '<');
				myFiles[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFiles[i].ignore(1024, ' ');
				myFiles[i] >> testStr;
			}
			
		}
		for (testStr = ""; testStr != "Name=\"id\"";  )
		{		 
			for (testStr = ""; testStr != "DataArray";  )
			{
				myFilesNext[i].ignore(1024, '<');
				myFilesNext[i] >> testStr;
			}
			for (testStr = ""; testStr.substr(0,4) != "Name";  )
			{
				myFilesNext[i].ignore(1024, ' ');
				myFilesNext[i] >> testStr;
			}
		}
		
		for (int j = 0; j < pointCount; j++)
		{
			myFiles[i].ignore(1024, '\n');
			myFilesNext[i].ignore(1024, '\n');
			myFiles[i] >> ids[i][j];
			myFilesNext[i] >> idsNext[i][j];
			cout << " id [" << i << "][" << j << 
			  "]: " << ids[i][j] <<
			  ", "  << idsNext[i][j] << "\n";
		}
		
		
		for (int j = 0; j < pointCount; j++)
		{
		  if (ids[i][j] == idsNext[i][j])
			  mapPointsNext[i][j] = j;
			else
			{
			  for (int k = 0; k < fpointCount; k++)
				{
				  if (ids[i][j] == idsNext[i][k])
					  mapPointsNext[i][j] = k;
				}
			}
		}
		
		
		for (int j = 0; (j < pointCount) && (j < fpointCount); j++)
		{
			cout << " Position [" << i << "][" << (remapPoints[i][j]) << 
			  "]: " << positions[i][(remapPoints[i][j])][0] <<
			  ", "  << positions[i][(remapPoints[i][j])][1] << 
			  ", "  << positions[i][(remapPoints[i][j])][2] <<
        ". Remap[" << i << "][" << j << "]: " << remapPoints[i][j] << "\n";
			cout << " Next Pos [" << i << "][" << (remapPoints[i][j]) << 
			  "]: " << positionsNext[i][(remapPoints[i][j])][0] <<
			  ", "  << positionsNext[i][(remapPoints[i][j])][1] << 
			  ", "  << positionsNext[i][(remapPoints[i][j])][2] <<
        ". Remap[" << i << "][" << j << "]: " << remapPoints[i][j] << "\n";
		}
		cout << "\n";
	}
	
	// end of reading in the points for each pair of files.o
	// the second file in the pair  has the final positions of the points.

	//Comparisons start here
	cout << "\n---Compare Begin---\n\n";

	for (int i = 0; i < fpointCount; i++) // loop over particles i and files j
	{

	  /*cout << "file[0]:Point[" << i << "]: <" << 
      positionsNext[0][mapPointsNext[0][i]][0] << "," <<
      positionsNext[0][mapPointsNext[0][i]][1] << "," <<
      positionsNext[0][mapPointsNext[0][i]][2] << "," <<

		  ">; ID: " << idsNext[0][mapPointsNext[0][i]] <<"\n";*/

		cout << "Next Point\n";

		int k2 = mapPointsNext[0][i]; // initially k2 is the index of the final position of the ith particle
		                              // in the first file pair, next time through the loop it is the
		                              // the final position of the ith particle	(j-1)st file pair

		double Last_Error = 0.0;

	  for (int j = 1; j < fileCount; j++)
		{

		  // set k1 to the start position to the ith particle

		  int k1 = remapPoints[j][i];

		  // set k1 to the index of the end position of the ith particle in the jth file pair

		  k1 = mapPointsNext[j][k1];

		    // int k2 = remapPoints[j-1][i];
			//     k2 = mapPointsNext[j-1][k1];


		  // Now calculate the error between the two particles

		  double errX = positionsNext[j][k1][0] - positionsNext[j-1][k2][0],
			       errY = positionsNext[j][k1][1] - positionsNext[j-1][k2][1],
			       errZ = positionsNext[j][k1][2] - positionsNext[j-1][k2][2];

		  double error = sqrt((errX * errX) + (errY * errY) + (errZ * errZ));

			
		  /*cout << " file[" << j-1 << "].point[" << k1 << "]: <" << 
				positionsNext[j-1][k2][0] << "," <<
				positionsNext[j-1][k2][1] << "," <<
				positionsNext[j-1][k2][2] << "," <<
				">; ID: " << idsNext[j-1][k2] << "\n";*/

		  // Now use Richardson Extrapolation to find the order of accuracy of this pair of particles
			
		 // Output this for one pair of particles, not multiple particles

			double Convergence_Rate = log2(error/Last_Error);
			cout << "                  Order: {" << Convergence_Rate << "}\n";
			
			cout << " [" << j-1 << "],[" << j << "]\n";
			 // ID: (" << idsNext[j-1][k2]
			     // << "," << idsNext[j][k1] << ")\n";
			cout << "  error: " << error << "\n";
			
			Last_Error = error;
			
		  /*cout << " file[" << j << "].point[" << k1 << "]: <" << 
				positionsNext[j][k1][0] << "," <<
				positionsNext[j][k1][1] << "," <<
				positionsNext[j][k1][2] << "," <<
				">; ID: " << idsNext[j][k1] << "\n";*/
			k2 = k1;
			//cout << "\n";
		}
		cout << "\n";
	}
	
	return 0;
}

///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <set>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face " << face / 3 << " ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << " " << otherHalf[dirEdge] << std::endl;
    } // WriteObjectStream()


//compute vertex's neighbor vertices, store them in neighborVertex
void DirectedEdgeSurface::computeNeighborVertex(unsigned int i) {
    neighborVertex[faceVertices[i]].push_back(faceVertices[i+1]);
    neighborVertex[faceVertices[i]].push_back(faceVertices[i+2]);
    neighborVertex[faceVertices[i+1]].push_back(faceVertices[i]);
    neighborVertex[faceVertices[i+1]].push_back(faceVertices[i+2]);
    neighborVertex[faceVertices[i+2]].push_back(faceVertices[i]);
    neighborVertex[faceVertices[i+2]].push_back(faceVertices[i+1]);
}
//update old positions
void DirectedEdgeSurface::updatePosition(unsigned int i) {
    vector<unsigned int> neighbours = neighborVertex[i];

    //delete the duplicated vertexID
    sort(neighbours.begin(), neighbours.end());
    auto it = unique(neighbours.begin(), neighbours.end());
    neighbours.erase(it, neighbours.end());

    //compute the sum of neighbor vertex
    float s = neighbours.size();
    float u = 1.0f / s * (5.0f / 8.0f - pow(3.0f / 8.0f + 1.0f / 4.0f * cos(2.0f * PI / s), 2));
    Cartesian3 pos_sum;
    for (unsigned int k = 0;k < s; k++)
    {
        pos_sum = pos_sum + vertices[neighbours[k]];
    }
    //update vertex position
    vertices.push_back(vertices[i] * (1.f -  s * u) + u * pos_sum);
}

//compute normals in a face using barycentre
void DirectedEdgeSurface::computeNormal(unsigned int i) {

    Cartesian3 v1 = vertices[faceVertices[i]];
    Cartesian3 v2 = vertices[faceVertices[i + 1]];
    Cartesian3 v3 = vertices[faceVertices[i + 2]];

    Cartesian3 v12 = v2 - v1;
    Cartesian3 v13 = v3 - v1;

    Cartesian3 v23 = v3 - v2;
    Cartesian3 v21 = v2 - v1;

    Cartesian3 v31 = v3 - v1;
    Cartesian3 v32 = v3 - v2;

    normals[faceVertices[i]] = v12.cross(v13).unit();
    normals[faceVertices[i + 1]] = v23.cross(v21).unit();
    normals[faceVertices[i + 2]] = v31.cross(v32).unit();
}

//compute next edge index in the sanme face
unsigned int DirectedEdgeSurface::next(unsigned int Edgeindex) {
    unsigned int face_index = Edgeindex / 3;
    unsigned int next = (Edgeindex % 3 + 1) % 3;
    return (face_index * 3 + next);
}

//main function
void DirectedEdgeSurface::loopSubdivision() {
    unsigned int vertex_num = vertices.size();
    unsigned int face_num = faceVertices.size();

        // loop every face
    for (unsigned int i = 0; i < face_num; i += 3){

            //compute and store neighbor vertex first
            computeNeighborVertex(i);

            //create a space to store new added vertex ID
            unsigned int new_vertexID[3];

            //store three vertices in a face and the opposite vertex
            unsigned int fourvertex[4];

            //first three vertex are vertices in a face
            for(unsigned int k=0;k<3;k++){
                fourvertex[k] = faceVertices[i+k];
            }

            for(unsigned int j = 0; j < 3; j++){

                //check if a new vertex is added at this edge or not
                if(!addFlag[i+j]){

                    //find opposite vertex in oppposite face
                    unsigned int other = otherHalf[i+j];
                    fourvertex[3] = faceVertices[next(other)];

                    //compute new position
                    Cartesian3 new_pos = (3.0f / 8.0f) * (vertices[fourvertex[j]] + vertices[fourvertex[next(next(j))]])
                             + (1.0f / 8.0f) * (vertices[fourvertex[next(j)]] + vertices[fourvertex[3]]);

                    //set the addFlag using edge index and vertexID
                    new_vertexID[j] = vertices.size();
                    addFlag[i+j] = new_vertexID[j];
                    addFlag[other] = new_vertexID[j];

                    //push new vertex and normal to origin vector
                    vertices.push_back((new_pos));
                    //normals.push_back(new_normal);
                }
                else{
                    //if a new vertex is added in this edge, find vertexID using edge index
                    new_vertexID[j] = addFlag[i+j];
                }
            }

            // add new top face
            faceVertices.push_back(faceVertices[i]);
            faceVertices.push_back(new_vertexID[1]);
            faceVertices.push_back(new_vertexID[0]);

            // add new left face
            faceVertices.push_back(faceVertices[i+1]);
            faceVertices.push_back(new_vertexID[2]);
            faceVertices.push_back(new_vertexID[1]);

            // add new right face
            faceVertices.push_back(faceVertices[i+2]);
            faceVertices.push_back(new_vertexID[0]);
            faceVertices.push_back(new_vertexID[2]);


            //update old face
            faceVertices[i] = new_vertexID[0];
            faceVertices[i+1] = new_vertexID[1];
            faceVertices[i+2] = new_vertexID[2];
        }

        unsigned int new_vertexnum = vertices.size();

        normals.clear();
        normals.resize(vertices.size());
        for(unsigned int i=0;i < faceVertices.size();i+=3){
            computeNormal(i);
        }
        // update the old vertices
        for (unsigned int i = 0; i < vertex_num; i++)
        {
            updatePosition(i);
        }
        unsigned int q = vertices.size();
        //put updated position to right place
        for(unsigned int i = new_vertexnum;i < q; i++){
            vertices[i - new_vertexnum] = vertices[i];
        }
        //delete the duplicate vertex
        vertices.erase(vertices.begin() + new_vertexnum ,vertices.end());

        //compute firstdirectedEdge
        firstDirectedEdge.clear();
        for(unsigned int i = 0; i< vertices.size(); i++){
            for(unsigned int j = 0; j < faceVertices.size(); j++){
                if(i == faceVertices[j]){
                    firstDirectedEdge.push_back(next(j));
                    break;
                }
            }
        }

        //compute other half
        otherHalf.clear();
        for(unsigned int i = 0; i < faceVertices.size(); i++){
            unsigned int vertexfrom = faceVertices[next(next(i))];
            unsigned int vertexto = faceVertices[i];
            for(unsigned int j = 0; j < faceVertices.size(); j++){
                if(vertexto == faceVertices[next(next(j))] && vertexfrom == faceVertices[j]){
                    otherHalf.push_back(j);
                }
            }
        }

        addFlag.clear();
        neighborVertex.clear();

}
// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting

    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()


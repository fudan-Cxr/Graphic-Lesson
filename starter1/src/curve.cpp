#include "curve.h"
#include "vertexrecorder.h"
using namespace std;

const float c_pi = 3.14159265358979323846f;

namespace
{
// Approximately equal to.  We don't want to use == because of
// precision issues with floating point.
inline bool approx(const Vector3f& lhs, const Vector3f& rhs)
{
	const float eps = 1e-8f;
	return (lhs - rhs).absSquared() < eps;
}


}


Curve evalBezier(const vector< Vector3f >& P, unsigned steps)
{
	// Check
	if (P.size() < 4 || P.size() % 3 != 1)
	{
		cerr << "evalBezier must be called with 3n+1 control points." << endl;
		exit(0);
	}

	// TODO:
	// You should implement this function so that it returns a Curve
	// (e.g., a vector< CurvePoint >).  The variable "steps" tells you
	// the number of points to generate on each piece of the spline.
	// At least, that's how the sample solution is implemented and how
	// the SWP files are written.  But you are free to interpret this
	// variable however you want, so long as you can control the
	// "resolution" of the discretized spline curve with it.

	// Make sure that this function computes all the appropriate
	// Vector3fs for each CurvePoint: V,T,N,B.
	// [NBT] should be unit and orthogonal.

	// Also note that you may assume that all Bezier curves that you
	// receive have G1 continuity.  Otherwise, the TNB will not be
	// be defined at points where this does not hold.


	


	// cerr << "\t>>> evalBezier has been called with the following input:" << endl;

	// cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
	// for (int i = 0; i < (int)P.size(); ++i)
	// {
	// 	cerr << "\t>>> " << P[i] << endl;
	// }

	// cerr << "\t>>> Steps (type steps): " << steps << endl;
	// cerr << "\t>>> Returning empty curve." << endl;

	// // Right now this will just return this empty curve.
	// return Curve();



	// Curve curve;
    
    // // 验证控制点数量 (3n+1)
    // if (P.size() < 4 || P.size() % 3 != 1)
    // {
    //     cerr << "evalBezier must be called with 3n+1 control points." << endl;
    //     exit(0);
    // }

    // // 计算每段贝塞尔曲线的控制点数量
    // const int segments = (P.size() - 1) / 3;
    
    // // 遍历每一段贝塞尔曲线
    // for (int i = 0; i < segments; ++i)
    // {
    //     // 获取当前段的4个控制点
    //     Vector3f p0 = P[3*i];
    //     Vector3f p1 = P[3*i + 1];
    //     Vector3f p2 = P[3*i + 2];
    //     Vector3f p3 = P[3*i + 3];
        
    //     // 为当前段生成曲线点
    //     for (unsigned j = 0; j <= steps; ++j)
    //     {
    //         float t = j / (float)steps;
            
    //         // 计算位置 (三次贝塞尔曲线公式)
    //         Vector3f position = 
    //             pow(1-t, 3) * p0 + 
    //             3 * pow(1-t, 2) * t * p1 + 
    //             3 * (1-t) * pow(t, 2) * p2 + 
    //             pow(t, 3) * p3;
            
    //         // 计算切线 (一阶导数)
    //         Vector3f tangent = 
    //             3 * pow(1-t, 2) * (p1 - p0) + 
    //             6 * (1-t) * t * (p2 - p1) + 
    //             3 * pow(t, 2) * (p3 - p2);
            
    //         // 归一化切线
    //         if (tangent.abs() > 0)
    //             tangent.normalize();
            
    //         // 计算二阶导数 (用于计算法线)
    //         Vector3f secondDerivative = 
    //             6 * (1-t) * (p2 - 2*p1 + p0) + 
    //             6 * t * (p3 - 2*p2 + p1);
            
    //         // 计算法线 (使用Frenet框架)
    //         Vector3f normal;
    //         if (secondDerivative.abs() > 0)
    //         {
    //             // 使用二阶导数减去切线方向的投影
    //             normal = secondDerivative - tangent * Vector3f::dot(tangent, secondDerivative);
    //             if (normal.abs() > 0)
    //                 normal.normalize();
    //         }
    //         else
    //         {
    //             // 如果二阶导数为零，选择任意垂直向量
    //             if (tangent[0] != 0 || tangent[1] != 0)
    //                 normal = Vector3f(-tangent[1], tangent[0], 0).normalized();
    //             else
    //                 normal = Vector3f(0, -tangent[2], tangent[1]).normalized();
    //         }
            
    //         // 计算副法线 (叉积)
    //         Vector3f binormal = Vector3f::cross(tangent, normal).normalized();
            
    //         // 确保法线垂直于切线
    //         normal = Vector3f::cross(binormal, tangent).normalized();
            
    //         // 创建曲线点
    //         CurvePoint point;
    //         point.V = position;
    //         point.T = tangent;
    //         point.N = normal;
    //         point.B = binormal;
            
    //         // 添加到曲线
    //         curve.push_back(point);
    //     }
    // }
    
    // return curve;



	int control_group = P.size() / 3;
	Curve Bez(control_group * steps);
	for (int i = 0; i < control_group; i++)
	{

		for (int j = 0; j < steps; j++)
		{
			float t = float(j) / steps;
			Bez[steps * i + j].V = (1 - t) * (1 - t) * (1 - t) * P[3 * i] + 3 * t * (1 - t) * (1 - t) * P[3 * i + 1] + 3 * t * t * (1 - t) * P[3 * i + 2] + t * t * t * P[3 * i + 3];

			Bez[steps * i + j].T = (-3 * (1 - t) * (1 - t) * P[3 * i] + 3 * (1 - 3 * t) * (1 - t) * P[3 * i + 1] + 3 * t * (2 - 3 * t) * P[3 * i + 2] + 3 * t * t * P[3 * i + 3]).normalized();

			if (i + j == 0)
				Bez[steps * i + j].N = Vector3f::cross(Vector3f(0, 0, 1), Bez[steps * i + j].T).normalized();
			else
				Bez[steps * i + j].N = Vector3f::cross(Bez[steps * i + j - 1].B, Bez[steps * i + j].T).normalized();

			Bez[steps * i + j].B = Vector3f::cross(Bez[steps * i + j].T, Bez[steps * i + j].N).normalized();
		}
	}

	// Right now this will just return this empty curve.
	return Bez;
}

Curve evalBspline(const vector< Vector3f >& P, unsigned steps)
{
	// Check
	if (P.size() < 4)
	{
		cerr << "evalBspline must be called with 4 or more control points." << endl;
		exit(0);
	}

	// TODO:
	// It is suggested that you implement this function by changing
	// basis from B-spline to Bezier.  That way, you can just call
	// your evalBezier function.

	int control_group = P.size() - 3;
	Curve Bsp(control_group * steps);

	for (int i = 0; i < control_group; i++)
	{

		for (int j = 0; j < steps; j++)
		{
			float t = float(j) / steps;
			Bsp[steps * i + j].V = (1 - t) * (1 - t) * (1 - t) / 6.0f * P[i] + (4 - 6 * t * t + 3 * t * t * t) / 6.0f * P[i + 1] + (1 + 3 * t + 3 * t * t - 3 * t * t * t) / 6.0f * P[i + 2] + (t * t * t) / 6.0f * P[i + 3];

			Bsp[steps * i + j].T = (-(1 - t) * (1 - t) / 2.0f * P[i] + (-4 * t + 3 * t * t) / 2.0f * P[i + 1] + (1 + 2 * t - 3 * t * t) / 2.0f * P[i + 2] + (t * t) / 2.0f * P[i + 3]).normalized();

			if (i + j == 0)
				Bsp[steps * i + j].N = Vector3f::cross(Vector3f(0, 0, 1), Bsp[steps * i + j].T).normalized();
			else
				Bsp[steps * i + j].N = Vector3f::cross(Bsp[steps * i + j - 1].B, Bsp[steps * i + j].T).normalized();

			Bsp[steps * i + j].B = Vector3f::cross(Bsp[steps * i + j].T, Bsp[steps * i + j].N).normalized();
		}
	}
	if (P[control_group] == P[0] && P[control_group + 1] == P[1] && P[control_group + 2] == P[2])
		Bsp.push_back(Bsp[0]);

	// Return an empty curve right now.
	return Bsp;




	// cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

	// cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
	// for (int i = 0; i < (int)P.size(); ++i)
	// {
	// 	cerr << "\t>>> " << P[i] << endl;
	// }

	// cerr << "\t>>> Steps (type steps): " << steps << endl;
	// cerr << "\t>>> Returning empty curve." << endl;

	// // Return an empty curve right now.
	// return Curve();
}

Curve evalCircle(float radius, unsigned steps)
{
	// This is a sample function on how to properly initialize a Curve
	// (which is a vector< CurvePoint >).

	// Preallocate a curve with steps+1 CurvePoints
	Curve R(steps + 1);

	// Fill it in counterclockwise
	for (unsigned i = 0; i <= steps; ++i)
	{
		// step from 0 to 2pi
		float t = 2.0f * c_pi * float(i) / steps;

		// Initialize position
		// We're pivoting counterclockwise around the y-axis
		R[i].V = radius * Vector3f(cos(t), sin(t), 0);

		// Tangent vector is first derivative
		R[i].T = Vector3f(-sin(t), cos(t), 0);

		// Normal vector is second derivative
		R[i].N = Vector3f(-cos(t), -sin(t), 0);

		// Finally, binormal is facing up.
		R[i].B = Vector3f(0, 0, 1);
	}

	return R;
}

void recordCurve(const Curve& curve, VertexRecorder* recorder)
{
	const Vector3f WHITE(1, 1, 1);
	for (int i = 0; i < (int)curve.size() - 1; ++i)
	{
		recorder->record_poscolor(curve[i].V, WHITE);
		recorder->record_poscolor(curve[i + 1].V, WHITE);
	}
}
void recordCurveFrames(const Curve& curve, VertexRecorder* recorder, float framesize)
{
	Matrix4f T;
	const Vector3f RED(1, 0, 0);
	const Vector3f GREEN(0, 1, 0);
	const Vector3f BLUE(0, 0, 1);
	
	const Vector4f ORGN(0, 0, 0, 1);
	const Vector4f AXISX(framesize, 0, 0, 1);
	const Vector4f AXISY(0, framesize, 0, 1);
	const Vector4f AXISZ(0, 0, framesize, 1);

	for (int i = 0; i < (int)curve.size(); ++i)
	{
		T.setCol(0, Vector4f(curve[i].N, 0));
		T.setCol(1, Vector4f(curve[i].B, 0));
		T.setCol(2, Vector4f(curve[i].T, 0));
		T.setCol(3, Vector4f(curve[i].V, 1));
 
		// Transform orthogonal frames into model space
		Vector4f MORGN  = T * ORGN;
		Vector4f MAXISX = T * AXISX;
		Vector4f MAXISY = T * AXISY;
		Vector4f MAXISZ = T * AXISZ;

		// Record in model space
		recorder->record_poscolor(MORGN.xyz(), RED);
		recorder->record_poscolor(MAXISX.xyz(), RED);

		recorder->record_poscolor(MORGN.xyz(), GREEN);
		recorder->record_poscolor(MAXISY.xyz(), GREEN);

		recorder->record_poscolor(MORGN.xyz(), BLUE);
		recorder->record_poscolor(MAXISZ.xyz(), BLUE);
	}
}


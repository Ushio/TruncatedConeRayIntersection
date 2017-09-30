#include "ofApp.h"

// Solve[a*sqrt((i + x*t)^2+(k + z*t)^2) + b*(j + y*t) + c == 0, t]

struct Cone {
	float h = 1.0f;
	float r = 1.0f;
	float min_h = 0.0f;
	float max_h = 1.0f;
};
inline float Sqr(float x) {
	return x * x;
}

inline bool intersect_cone(glm::vec3 o, glm::vec3 d, Cone cone, float *tmin) {
	float h = cone.h;
	float r = cone.r;

	float k = Sqr(h) / Sqr(r);

	float a = k * (d.x * d.x + d.z * d.z) - d.y * d.y;
	float b = 2.0f * k * (o.x * d.x + o.z * d.z) - 2.0f * (o.y - h) * d.y;
	float c = k * (o.x * o.x + o.z * o.z) - Sqr(o.y - h);

	float D = b * b - 4.0f * a * c;
	if (D < 0.0f) {
		return false;
	}

	float t;
	bool isect = false;

	t = (-b + sqrt(D)) / (2.0f * a);
	if (0.0f < t && t < *tmin) {
		float y = o.y + d.y * t;
		if (cone.min_h <= y && y <= cone.max_h) {
			*tmin = t;
			isect = true;
		}
	}

	t = (-b - sqrt(D)) / (2.0f * a);
	if (0.0f < t && t < *tmin) {
		float y = o.y + d.y * t;
		if (cone.min_h <= y && y <= cone.max_h) {
			*tmin = t;
			isect = true;
		}
	}

	return isect;
}

// 円錐の方程式の、(∂u/∂x, ∂u/∂y, ∂u/∂z)
inline glm::vec3 gradient_cone(glm::vec3 p, Cone cone) {
	// D[(h*x/r)^2 +(h*z/r)^2 - (y-h)^2, x]
	float dudx = 2.0f * cone.h * cone.h * p.x / (cone.r * cone.r);
	float dudy = 2.0f * (cone.h - p.y);
	float dudz = 2.0f * cone.h * cone.h * p.z / (cone.r * cone.r);
	return glm::vec3(dudx, dudy, dudz);
}

struct Cylinder {
	float r = 1.0f;
	float min_h = 0.0f;
	float max_h = 1.0f;
};
inline bool intersect_cylinder(glm::vec3 o, glm::vec3 d, Cylinder cylinder, float *tmin) {
	float r = cylinder.r;

	float a = d.x * d.x + d.z * d.z;
	float b = 2.0f * (o.x * d.x + o.z * d.z);
	float c = o.x * o.x + o.z * o.z - r * r;

	float D = b * b - 4.0f * a * c;
	if (D < 0.0f) {
		return false;
	}

	float t;
	bool isect = false;

	t = (-b + sqrt(D)) / (2.0f * a);
	if (0.0f < t && t < *tmin) {
		float y = o.y + d.y * t;
		if (cylinder.min_h <= y && y <= cylinder.max_h) {
			*tmin = t;
			isect = true;
		}
	}

	t = (-b - sqrt(D)) / (2.0f * a);
	if (0.0f < t && t < *tmin) {
		float y = o.y + d.y * t;
		if (cylinder.min_h <= y && y <= cylinder.max_h) {
			*tmin = t;
			isect = true;
		}
	}

	return isect;
}


struct TruncatedCone {
	glm::vec3 p;
	float p_radius = 0.0f;
	glm::vec3 q;
	float q_radius = 1.0f;
};
inline bool intersect_cone(glm::vec3 o, glm::vec3 d, TruncatedCone cone, float *tmin, glm::vec3 *n) {
	float L = glm::distance(cone.p, cone.q);
	if (L <= glm::epsilon<float>()) {
		return false;
	}

	if (cone.q_radius < cone.p_radius) {
		std::swap(cone.p, cone.q);
		std::swap(cone.p_radius, cone.q_radius);
	}
	
	// p is upper
	glm::vec3 cone_o = cone.q;
	glm::vec3 cone_d = (cone.p - cone.q) / L;
	
	// 円錐を水平に戻す
	// ここはライブラリによってはcone_d.y < -0.999...のとき、未定義になる可能性あり
	glm::quat r_to_cone_local = glm::rotation(cone_d, glm::vec3(0.0f, 1.0f, 0.0f));
	glm::vec3 t_to_cone_local = -cone_o;
	glm::vec3 d_cone_local = r_to_cone_local * d;
	glm::vec3 o_cone_local = r_to_cone_local * (o + t_to_cone_local);

	if (std::fabs(cone.q_radius - cone.p_radius) < 0.000001f) {
		// intersect as cylinder
		Cylinder cylinder;
		cylinder.r = (cone.q_radius + cone.p_radius) * 0.5f;
		cylinder.min_h = 0.0f;
		cylinder.max_h = L;
		if (intersect_cylinder(o_cone_local, d_cone_local, cylinder, tmin)) {
			glm::vec3 p = o_cone_local + d_cone_local * (*tmin);
			*n = glm::normalize(glm::vec3(p.x, 0.0f, p.z));
			if (0.0f < glm::dot(*n, d_cone_local)) {
				*n = -*n;
			}
			*n = glm::inverse(r_to_cone_local) * (*n);
			return true;
		} else {
			return false;
		}
	}

	// intersect as cone

	Cone pure_cone;
	pure_cone.h = cone.q_radius / (cone.q_radius - cone.p_radius) * L;
	pure_cone.r = cone.q_radius;
	pure_cone.min_h = 0.0f;
	pure_cone.max_h = L;
	if (intersect_cone(o_cone_local, d_cone_local, pure_cone, tmin)) {
		glm::vec3 p = o_cone_local + d_cone_local * (*tmin);

		/*
		三角錐の方程式 ~~~~ = 0 において、
		(∂u/∂x)⊿x + (∂u/∂y)⊿y + (∂u/∂z)⊿z == 0
		となるのは、(⊿x, ⊿y, ⊿z) の移動が三角錐の表面上である場合のみである。
		その移動ベクトルは面をなし、その移動ベクトルと(∂u/∂x, ∂u/∂y, ∂u/∂z)の内積がゼロである。
		これはつまり、(∂u/∂x, ∂u/∂y, ∂u/∂z)は表面と直交している。よって、これは法線ベクトルとなる。
		*/
		*n = glm::normalize(gradient_cone(p, pure_cone));
		if (0.0f < glm::dot(*n, d_cone_local)) {
			*n = -*n;
		}
		*n = glm::inverse(r_to_cone_local) * (*n);
		return true;
	}
	return false;
}

inline ofVec3f toOf(glm::vec3 p) {
	return ofVec3f(p.x, p.y, p.z);
}

//--------------------------------------------------------------
void ofApp::setup(){
	ofSetVerticalSync(false);
	ofSetFrameRate(60.0f);

	_imgui.setup();

	_camera.setNearClip(0.1f);
	_camera.setFarClip(100.0f);
	_camera.setDistance(5.0f);

	ofSeedRandom(0);
}

//--------------------------------------------------------------
void ofApp::update() {

}

inline void drawCone(TruncatedCone cone) {
	if (cone.q_radius < cone.p_radius) {
		std::swap(cone.p, cone.q);
		std::swap(cone.p_radius, cone.q_radius);
	}
	// p is upper
	glm::vec3 cone_o = cone.q;
	glm::vec3 cone_d = glm::normalize(cone.p - cone.q);

	// ここはライブラリによってはcone_d.y < -0.999...のとき、未定義になる可能性あり
	glm::quat r_to_cone_local = glm::rotation(cone_d, glm::vec3(0.0f, 1.0f, 0.0f));
	glm::vec3 t_to_cone_local = -cone_o;
	glm::mat4 m_to_cone_local = glm::mat4_cast(r_to_cone_local) * glm::translate(t_to_cone_local);

	ofPushStyle();
	ofNoFill();

	ofPushMatrix();
	ofMultMatrix(glm::value_ptr(glm::inverse(m_to_cone_local)));

	float L = glm::distance(cone.p, cone.q);

	int N = 10;
	for (int i = 0; i < N; ++i) {
		float theta = glm::two_pi<float>() * i / float(N);
		float rx = cos(theta);
		float rz = sin(theta);
		ofLine(
			cone.q_radius * rx, 0, cone.q_radius * rz,
			cone.p_radius * rx, L, cone.p_radius * rz
		);
	}

	ofPushMatrix();
	ofTranslate(0, L, 0);
	ofRotateX(90);
	ofDrawCircle(0, 0, cone.p_radius);
	ofPopMatrix();

	ofPushMatrix();
	ofRotateX(90);
	ofDrawCircle(0, 0, cone.q_radius);
	ofPopMatrix();

	ofPopMatrix();

	ofPopStyle();
}

//--------------------------------------------------------------
void ofApp::draw(){

	// drawing
	ofEnableDepthTest();

	ofClear(0);

	_camera.begin();

	// floor
	ofPushMatrix();
	ofRotateZ(90.0f);
	ofSetColor(64);
	ofDrawGridPlane(1.0f);
	ofPopMatrix();


	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	static TruncatedCone cone;
	static std::once_flag cone_init_flag;

	std::call_once(cone_init_flag, []() {
		cone.p = glm::vec3(0.0f, 0.0f, 1.0f);
		cone.p_radius = 0.2f;
		cone.q = glm::vec3(0.0f, 1.0f, 1.0f);
		cone.q_radius = 1.0f;
	});

	ofSetColor(255);
	drawCone(cone);

	static float dx = 0.0f;
	static float Z = 3.0f;

	static int focus_x = 0;
	static int focus_y = 0;

	int N = 50;
	for (int i = 0; i < N; ++i) {
		float s0 = (float)i / (N - 1);
		float x = ofLerp(-1.5f, 1.5f, s0);

		for (int j = 0; j < N; ++j) {
			float s1 = (float)j / (N - 1);
			float y = ofLerp(-1.5f, 1.5f, s1);

			if (focus_x == i && focus_y == j) {
				printf("");
			}

			glm::vec3 orig(x, y, Z);
			glm::vec3 dir = glm::normalize(glm::vec3(dx, 0.0f, -1.0f));

			float tmin = 1000000.0f;
			glm::vec3 n;
			if (intersect_cone(orig, dir, cone, &tmin, &n)) {
				auto p = orig + dir * tmin;

				ofSetColor(255, 125, 0);
				ofLine(toOf(orig), toOf(p));
				ofDrawSphere(toOf(p), 0.01f);

				ofSetColor(64, 125, 255);
				// ofLine(toOf(p), toOf(p + glm::reflect(dir, n) * 0.2f));
				ofLine(toOf(p), toOf(p + n * 0.2f));
			}
			else {
				if (focus_x == i && focus_y == j) {
					ofSetColor(255, 64, 64, 128);
				}
				else {
					ofSetColor(64, 64, 64);
				}
				ofLine(toOf(orig), toOf(orig + dir * 50.0f));
			}
		}
	}

	_camera.end();


	ofDisableDepthTest();
	ofSetColor(255);

	_imgui.begin();

	ImGui::PushStyleColor(ImGuiCol_WindowBg, ofVec4f(0.2f, 0.2f, 0.5f, 0.5f));
	ImGui::SetNextWindowPos(ofVec2f(500, 30), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ofVec2f(500, 600), ImGuiSetCond_Once);

	ImGui::Begin("Config Panel");
	ImGui::Text("fps: %.2f", ofGetFrameRate());

	ImGui::SliderFloat("dx", &dx, -1.0f, 1.0f);
	ImGui::SliderFloat("Z", &Z, -3.0f, 3.0f);
	ImGui::InputInt("focus x", &focus_x);
	ImGui::InputInt("focus y", &focus_y);
	
	if (ImGui::Button("random")) {
		cone.p = glm::vec3(ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f));
		cone.p_radius = ofRandom(0.0f, 2.0f);
		cone.q = glm::vec3(ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f));
		cone.q_radius = ofRandom(0.0f, 2.0f);
	}
	if (ImGui::Button("random as cylinder")) {
		cone.p = glm::vec3(ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f));
		cone.p_radius = ofRandom(0.0f, 2.0f);
		cone.q = glm::vec3(ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f), ofRandom(-1.0f, 1.0f));
		cone.q_radius = cone.p_radius;
	}

	auto wp = ImGui::GetWindowPos();
	auto ws = ImGui::GetWindowSize();
	ofRectangle win(wp.x, wp.y, ws.x, ws.y);

	ImGui::End();
	ImGui::PopStyleColor();

	_imgui.end();

	if (win.inside(ofGetMouseX(), ofGetMouseY())) {
		_camera.disableMouseInput();
	}
	else {
		_camera.enableMouseInput();
	}
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}

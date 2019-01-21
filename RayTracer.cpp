// https://habr.com/ru/post/436790/
// https://github.com/ssloy/tinyraytracer

#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "Geometry.h"

struct Light {
public:
    Light(const Vec3f &p, const float &i): 
		position(p), intensity(i) {
	}
public:
    Vec3f position;
    float intensity;
};

struct Material {
public:
    Material(const float &r, const Vec4f &a, const Vec3f &color, const float &spec): 
		refractive_index(r), 
		albedo(a), 
		diffuse_color(color), 
		specular_exponent(spec) {
	}
    Material(): 
		refractive_index(1), 
		albedo(1,0,0,0), 
		diffuse_color(), 
		specular_exponent() {
	}
public:
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

struct Sphere {
public:
    Vec3f center;
    float radius;
    Material material;

public:
    Sphere(const Vec3f &c, const float &r, const Material &m): 
		center(c), 
		radius(r), 
		material(m) {
	}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

Vec3f reflect(const Vec3f& I, const Vec3f& N) {
    return I - N*2.0f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index) { // Snell's law
    float cosi = -std::max(-1.f, std::min(1.f, I*N));
    float etai = 1, etat = refractive_index;
    Vec3f n = N;
    if (cosi < 0) { // if the ray is inside the object, swap the indices and invert the normal to get the correct result
        cosi = -cosi;
        std::swap(etai, etat); n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));
}

// Проверка пересечения со сферой
bool sceneIntersect(const Vec3f& orig, // Откуда мы бросаем луч
                    const Vec3f& dir,  // Направление
                    const std::vector<Sphere>& spheres, // Список сфер
                    Vec3f& hit, // Место, где произошло пересечение со сферой
                    Vec3f& N,   // Нормаль
                    Material& material) { // Материал сферы с которой произошло пересечение
    // Дистанция до ближайшей сцены
    float spheresMaxDist = std::numeric_limits<float>::max();
    for (size_t i = 0; i < spheres.size(); i++) {
        // Проверяем, есть ли пересечение со сферой
        float distI = 0; // Дистанция до сферы
        bool hasIntersect = spheres[i].ray_intersect(orig, dir, distI);
        
        // Было ли пересейчение со сферой и является ли сфера более близкой, чем предыдущая
        if (hasIntersect && (distI < spheresMaxDist)) {
            spheresMaxDist = distI; // Сохраняем дистанцию до сферы как максимальную
            hit = orig + dir*distI; // Точка пересечения - точка из экрана + направление, умноженное на расстояние до пересечения со сферой
            N = (hit - spheres[i].center).normalize(); // Нормаль в точке пересечения
            material = spheres[i].material; // Материал
        }
    }
    
    // TODO: ???
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 1e-3)  {
        float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheresMaxDist) {
            checkerboardDist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(1,1,1) : Vec3f(1, .7, .3);
            material.diffuse_color = material.diffuse_color*.3;
        }
    }
    return std::min(spheresMaxDist, checkerboardDist)<1000;
}

Vec3f castRay(const Vec3f& orig, // Откуда бросаем луч
              const Vec3f &dir,  // Направление луча
              const std::vector<Sphere>& spheres,   // Список сфер
              const std::vector<Light>& lights,     // Источники света
              size_t depth=0) { // Текущая глубина рекурсии
    
    // Ограничение глубины рекурсии
	const size_t depthLimit = 4;
	
    Vec3f point; // Точка пересечения со сферой
	Vec3f N;     // Нормаль в точке пересечения
    Material material; // Материал в точке пересечения

    if ((depth > depthLimit) || !sceneIntersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8); // background color
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = castRay(reflect_orig, reflect_dir, spheres, lights, depth + 1);
    Vec3f refract_color = castRay(refract_orig, refract_dir, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (sceneIntersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

void render(const std::vector<Sphere>& spheres, const std::vector<Light> &lights) {
	// Создаем буффер нужного размера под изображение
    const int width = 1024;
    const int height = 768;
    const int fov = M_PI/2.0; // Угол обзора камеры
    std::vector<Vec3f> framebuffer(width*height);
	
	const float imageRatio = width / (float)height;
	const float tanValue = tan(fov/2.0);
    
    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
			// Вычисляем нормализованное направление, по которому мы должны бросать луч из картинци в сцену
            float x = (2*(i + 0.5)/(float)width - 1) * tanValue * imageRatio;
            float y = -(2*(j + 0.5)/(float)height - 1) * tanValue;
            Vec3f dir = Vec3f(x, y, -1).normalize();
			
			// Пускаем луч
			Vec3f colorValue = castRay(Vec3f(0,0,0), dir, spheres, lights);
			
			// Сохраняем значение цвета
            framebuffer[i+j*width] = colorValue;
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
	// Создаем материалы
    Material ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), 50.0);
    Material glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), 125.0);
    Material redRubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1), 10.0);
    Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.0);
	
	// Создаем сферы
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, glass));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, redRubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4, mirror));
	
	// Создаем источники света
    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));
	
	// Выполняем рендеринг
    render(spheres, lights);

    return 0;
}


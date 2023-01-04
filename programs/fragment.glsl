#version 330 core
#include hg_sdf.glsl
layout (location = 0) out vec4 fragColor;


in vec2 uv_0;
uniform vec2 u_resolution;
uniform sampler2D u_texture_0;
uniform float count;
uniform vec2 u_mouse;

const int MAXDEPTH = 100;

float random(float x)
{
    float y = fract(sin(x)*100000.0);
    return y;
}
float rand(float co){
    return fract(sin( (co)*12.9898 * 43758.5453));
}

vec3 random(vec3 v){
    return vec3(random(v.x+count*7.21),random(v.y+count*24.341),random(v.z*count*2.546156));
}

vec3 random_in_unit_sphere(vec3 v){
    
    vec3 p =random(v);
        
    return normalize(p);
    

}



mat3 getCam(vec3 ro,vec3 lookAt){
    vec3 camF = normalize(lookAt-ro);
    vec3 camR = normalize(cross(vec3(0,1,0),camF));
    vec3 camU = cross(camF,camR);
    return mat3(camR,camU,camF);

}
void mouseControll(inout vec3 ro){
    vec2 m = u_mouse/u_resolution;
    pR(ro.yz,m.y * PI * 0.5-0.5);
    pR(ro.xz,m.x*TAU);

}

struct Ray{
    vec3 oringin;
    vec3 dir;
    vec3 at(float t){
        return oringin + t*dir;
    }
};




struct material{
    int type;
    vec3 ambient;
    float rr;
    vec3 get_color(){
        if(type==2) return vec3(1.0,1.0,1.0);
        return ambient;
    }

    Ray scatter(vec3 p,vec3 normal,vec3 dir){
        Ray ray;
        if(type==0){    
            vec3 target = p + normal + random_in_unit_sphere(p);
            ray.oringin = p;
            ray.dir = target - p;
        }else if(type==1){
            vec3 r = reflect(normalize(dir),normal);
            ray.oringin = p;
            ray.dir = r;
        }
        return ray;
    }

};

struct hit_record{
    vec3 p;
    vec3 normal;
    float t;
    material m;
    bool front_face;

    void set_face_normal(Ray r,vec3 n){
        front_face = dot(r.dir,n) < 0;
        normal = front_face? n : -n;

    }
};


struct sphere{
    vec3 center;
    float radius;
    material m;
    bool hit(Ray r,float t_min,float t_max,inout hit_record rec){   
        vec3 oc = r.oringin - center;
        float a = dot(r.dir, r.dir);
        float b = 2.0 * dot(oc, r.dir);
        float c = dot(oc, oc) - radius*radius;
        float discriminant = b*b - 4*a*c;
        
        if(discriminant<0){
            return false;
        }
        float root = (-b - sqrt(discriminant))/2*a;
        if (root < t_min || t_max < root) {
            root = (-b + sqrt(discriminant)) / 2*a;
            if (root < t_min || t_max < root)
                return false;
        }

        
        rec.t = root;
        rec.p = r.at(root);
        vec3 normal = (rec.p - center)/radius;
        
        rec.set_face_normal(r, normal);
        rec.m = m;
        return true;

    }
};

sphere spheres[445];

vec3 myrefract(vec3 uv, vec3 n, float etai_over_etat) {
    float cos_theta = min(dot(-uv, n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -sqrt(abs(1.0 - dot(r_out_perp,r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}
float reflectance(float cosine, float ref_idx) {
    // Use Schlick's approximation for reflectance.
    float r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}

vec3 ray_color(Ray ray){
    
    hit_record rec;
    bool hit_anything = false;
    float csf = 1000000;
    float t = (ray.dir.y+1)*0.5;
    vec3 temp_col = (1-t)*vec3(1.0,1.0,1.0)+t*vec3(0.5,0.7,1.0);
    Ray temp_ray = ray;
    float d;

    for(int k=0;k<MAXDEPTH;k+=1){
        bool hit_anything = false;
        float csf = 100000;
        if(k==MAXDEPTH-1){
            temp_col*=vec3(0,0,0);
            break;
        }

        for(int i=0;i<spheres.length;i+=1){
            if(spheres[i].hit(ray,0.001,csf,rec)){
                hit_anything = true;        
                csf = rec.t;
            }

        }
        if(hit_anything){
            //ray = rec.m.scatter(rec.p,rec.normal,ray.dir);
            if(rec.m.type==1){
                vec3 r = reflect(normalize(ray.dir),rec.normal);
                ray.oringin = rec.p;
                ray.dir = r;
            }else if(rec.m.type==0){
                vec3 target = rec.p + rec.normal + random_in_unit_sphere(rec.p);
                ray.oringin = rec.p;
                ray.dir = target - rec.p;

            }else{
                float rr = rec.front_face?(1.0/rec.m.rr):rec.m.rr;
                vec3 unit_direction = normalize(ray.dir);
                float cos_theta = min(dot(-unit_direction, rec.normal), 1.0);
                float sin_theta = sqrt(1.0 - cos_theta*cos_theta);
                bool cannot_refract = rr * sin_theta > 1.0;
                vec3 direction;
                if (cannot_refract||reflectance(cos_theta,rr)>(random(count)+1)*0.5){
                    direction = reflect(unit_direction, rec.normal);
                }else{
                    direction = myrefract(unit_direction, rec.normal, rr);
                }
                
                ray.oringin = rec.p;
                ray.dir = direction;
            }
            temp_col*=rec.m.get_color()*vec3(0.6,0.6,0.6);      
        }else{
            break;
        }
    }
    
    vec3 col = temp_col;
    return col;
}



void render(inout vec3 col,in vec2 uv){
    vec3 ro = vec3(13,2,3);
    mouseControll(ro);
    vec3 lookAt = vec3(0,0,0);
    vec3 rd = getCam(ro,lookAt)*normalize(vec3(uv,1));
    vec3 rl = normalize(vec3(uv,-1)-ro);
    Ray r = Ray(ro,rd);
    col = ray_color(r);
    
    
}

void set_world(){
    
    sphere s1 = sphere(vec3(0,-1000,0),1000,material(0,vec3(0.5,0.5,0.5),0));    
    sphere s2 = sphere(vec3(0,1,0),1.0,material(2,vec3(0.8,0.8,0.0),1.5));
    sphere s3 = sphere(vec3(-4,1,0),1.0,material(0,vec3(0.4,0.2,0.1),0));
    sphere s4 = sphere(vec3(4,1,0),1.0,material(1,vec3(0.7,0.6,0.5),0));
    sphere s5 = sphere(vec3(-1,0,-1),-0.4,material(2,vec3(0.8,0.6,0.2),1.5));
    int c = 4;
    for(int i=-10;i<11;i+=1){
        for(int j=-10;j<11;j+=1){
            vec3 center = vec3(i,0.35,j);
            vec3 albedo = vec3(1.0,1.0,1.0);//random(vec3(i*5,j*2.3,i*j*2.13));
            spheres[c] = sphere(center,0.35,material(1,albedo,0));
            c+=1;
        }

    }
    
    spheres[0] = s1;
    spheres[1] = s2;
    spheres[2] = s3;
    spheres[3] = s4;
    
   
}

void main() {
    vec2 uv = (2.0 * gl_FragCoord.xy - u_resolution.xy)/u_resolution.y;
    vec3 col;
    vec3 color = texture(u_texture_0,uv_0,0).rgb ;

    set_world();
    uv += (vec2(random(1+count*3),random(100+count*5))+1)*0.5/u_resolution.y;
    render(col, uv);
    col = pow(col,vec3(0.4545));
    col = (col + color*count)/(count+1);
    
    fragColor = vec4(col, 1.0);
}
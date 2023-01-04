import moderngl_window as mglw
from moderngl_window import screenshot
import numpy as np


class App(mglw.WindowConfig):
    window_size = 1600,900 
    resource_dir = 'programs'
    count = 0

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        self.quad = mglw.geometry.quad_fs()
        self.texture = self.ctx.texture(self.window_size,components=3,data=self.ctx.screen.read(components=3))
        self.program = self.load_program(vertex_shader='vertex.glsl', fragment_shader='fragment.glsl')
        # uniforms
        self.program['u_resolution'] = self.window_size
        self.program['u_texture_0'] = 0
        self.texture.use()
        self.program['count'] = self.count
        

    def render(self, time, frame_time):
        self.ctx.clear()          
        self.quad.render(self.program)       
        self.texture = self.ctx.texture(self.window_size,components=3,data=self.ctx.screen.read(components=3))
     
        self.count+=1
        self.program['u_texture_0'] = 0
        self.texture.use()
        self.program['count'] = self.count
        
        
        
        

    def mouse_position_event(self, x: int, y: int, dx: int, dy: int):
        self.program['u_mouse'] = (x,y)
        self.count = 0


if __name__ == '__main__':
    mglw.run_window_config(App)
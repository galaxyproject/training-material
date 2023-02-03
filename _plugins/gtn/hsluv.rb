# The MIT License (MIT)
# 
# Copyright (c) 2015 Radu-Bogdan Croitoru
# Copyright (c) 2016 Alexei Boronine
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# https://github.com/hsluv/hsluv-ruby


module Hsluv
  extend self

  M = [
      [3.240969941904521, -1.537383177570093, -0.498610760293],
      [-0.96924363628087, 1.87596750150772, 0.041555057407175],
      [0.055630079696993, -0.20397695888897, 1.056971514242878],
  ]

  M_INV = [
      [0.41239079926595, 0.35758433938387, 0.18048078840183],
      [0.21263900587151, 0.71516867876775, 0.072192315360733],
      [0.019330818715591, 0.11919477979462, 0.95053215224966],
  ]

  REF_X = 0.95045592705167
  REF_Y = 1.0
  REF_Z = 1.089057750759878
  REF_U = 0.19783000664283
  REF_V = 0.46831999493879
  KAPPA = 903.2962962
  EPSILON = 0.0088564516

  ###

  def hsluv_to_hex (h, s, l)
    rgb_to_hex(*hsluv_to_rgb(h, s, l))
  end

  def hpluv_to_hex (h, s, l)
    rgb_to_hex(*hpluv_to_rgb(h, s, l))
  end

  def hex_to_hsluv (hex)
    rgb_to_hsluv(*hex_to_rgb(hex))
  end

  def hex_to_hpluv (hex)
    rgb_to_hpluv(*hex_to_rgb(hex))
  end

  def hsluv_to_rgb (h, s, l)
    xyz_to_rgb(luv_to_xyz(lch_to_luv(hsluv_to_lch([h, s, l]))))
  end

  def rgb_to_hsluv (r, g, b)
    lch_to_hsluv(rgb_to_lch(r, g, b))
  end

  def hpluv_to_rgb (h, s, l)
    lch_to_rgb(*hpluv_to_lch([h, s, l]))
  end

  def rgb_to_hpluv (r, g, b)
    lch_to_hpluv(rgb_to_lch(r, g, b))
  end

  def lch_to_rgb (l, c, h)
    xyz_to_rgb(luv_to_xyz(lch_to_luv([l, c, h])))
  end

  def rgb_to_lch (r, g, b)
    luv_to_lch(xyz_to_luv(rgb_to_xyz([r, g, b])))
  end

  def rgb_to_hex (r, g, b)
    '#%02x%02x%02x' % rgb_prepare([r, g, b])
  end

  def hex_to_rgb (hex)
    hex = hex.tr('#', '')
    [].tap { |arr| hex.split('').each_slice(2) { |block| arr << block.join.to_i(16) / 255.0 } }
  end

  ###

  def rgb_to_xyz (arr)
    rgbl = arr.map { |val| to_linear(val) }
    M_INV.map { |i| dot_product(i, rgbl) }
  end

  def xyz_to_luv (arr)
    x, y, z = arr
    l = f(y)

    return [0.0, 0.0, 0.0] if [x, y, z, 0.0].uniq.length == 1 || l == 0.0

    var_u = (4.0 * x) / (x + (15.0 * y) + (3.0 * z))
    var_v = (9.0 * y) / (x + (15.0 * y) + (3.0 * z))
    u = 13.0 * l * (var_u - REF_U)
    v = 13.0 * l * (var_v - REF_V)

    [l, u, v]
  end

  def luv_to_lch (arr)
    l, u, v = arr
    c = ((u ** 2) + (v ** 2)) ** (1 / 2.0)
    hrad = Math.atan2(v, u)
    h = radians_to_degrees(hrad)
    h += 360.0 if h < 0.0
    [l, c, h]
  end

  def lch_to_hsluv (arr)
    l, c, h = arr
    return [h, 0.0, 100.0] if l > 99.9999999
    return [h, 0.0, 0.0] if l < 0.00000001

    mx = max_chroma_for(l, h)
    s = c / mx * 100.0

    [h, s, l]
  end

  def lch_to_hpluv (arr)
    l, c, h = arr

    return [h, 0.0, 100.0] if l > 99.9999999
    return [h, 0.0, 0.0] if l < 0.00000001

    mx = max_safe_chroma_for(l)
    s = c / mx * 100.0

    [h, s, l]
  end

  ###

  def xyz_to_rgb (arr)
    xyz = M.map { |i| dot_product(i, arr) }
    xyz.map { |i| from_linear(i) }
  end

  def luv_to_xyz (arr)
    l, u, v = arr

    return [0.0, 0.0, 0.0] if l == 0

    var_y = f_inv(l)
    var_u = u / (13.0 * l) + REF_U
    var_v = v / (13.0 * l) + REF_V


    y = var_y * REF_Y
    x = 0.0 - (9.0 * y * var_u) / ((var_u - 4.0) * var_v - var_u * var_v)
    z = (9.0 * y - (15.0 * var_v * y) - (var_v * x)) / (3.0 * var_v)

    [x, y, z]
  end

  def lch_to_luv (arr)
    l, c, h = arr

    hrad = degrees_to_radians(h)
    u = Math.cos(hrad) * c
    v = Math.sin(hrad) * c

    [l, u, v]
  end

  def hsluv_to_lch (arr)
    h, s, l = arr

    return [100, 0.0, h] if l > 99.9999999
    return [0.0, 0.0, h] if l < 0.00000001

    mx = max_chroma_for(l, h)
    c = mx / 100.0 * s

    [l, c, h]
  end

  def hpluv_to_lch (arr)
    h, s, l = arr

    return [100, 0.0, h] if l > 99.9999999
    return [0.0, 0.0, h] if l < 0.00000001

    mx = max_safe_chroma_for(l)
    c = mx / 100.0 * s

    [l, c, h]
  end

  ###

  def radians_to_degrees (rad)
    rad * 180.0 / Math::PI
  end

  def degrees_to_radians (degrees)
    degrees * Math::PI / 180.0
  end

  def max_chroma_for (l, h)
    hrad = h / 360.0 * Math::PI * 2.0
    lengths = []

    get_bounds(l).each do |line|
      l = length_of_ray_until_intersect(hrad, line)
      lengths << l if l
    end

    lengths.min
  end

  def max_safe_chroma_for (l)
    lengths = []

    get_bounds(l).each do |m1, b1|
      x = intersect_line_line([m1, b1], [-1.0 / m1, 0.0])
      lengths << distance_from_pole([x, b1 + x * m1])
    end

    lengths.min
  end

  def get_bounds (l)
    sub1 = ((l + 16.0) ** 3.0) / 1560896.0
    sub2 = sub1 > EPSILON ? sub1 : l / KAPPA
    ret = []

    M.each do |m1, m2, m3|
      [0, 1].each do |t|
        top1 = (284517.0 * m1 - 94839.0 * m3) * sub2
        top2 = (838422.0 * m3 + 769860.0 * m2 + 731718.0 * m1) * l * sub2 - 769860.0 * t * l
        bottom = (632260.0 * m3 - 126452.0 * m2) * sub2 + 126452.0 * t
        ret << [top1 / bottom, top2 / bottom]
      end
    end

    ret
  end

  def length_of_ray_until_intersect (theta, line)
    m1, b1 = line
    length = b1 / (Math.sin(theta) - m1 * Math.cos(theta))
    return nil if length < 0
    length
  end

  def intersect_line_line (line1, line2)
    (line1[1] - line2[1]) / (line2[0] - line1[0])
  end

  def distance_from_pole (point)
    Math.sqrt(point[0] ** 2 + point[1] ** 2)
  end

  def f (t)
    t > EPSILON ? 116 * ((t / REF_Y) ** (1.0 / 3.0)) - 16.0 : t / REF_Y * KAPPA
  end

  def f_inv (t)
    t > 8 ? REF_Y * ((t + 16.0) / 116.0) ** 3.0 : REF_Y * t / KAPPA
  end

  def to_linear (c)
    c > 0.04045 ? ((c + 0.055) / 1.055) ** 2.4 : c / 12.92
  end

  def from_linear (c)
    c <= 0.0031308 ? 12.92 * c : (1.055 * (c ** (1.0 / 2.4)) - 0.055)
  end

  def dot_product (a, b)
    a.zip(b).map { |i, j| i * j }.inject(:+)
  end

  def rgb_prepare (arr)
    arr.map! { |ch| ch = ch.round(3); ch = [0, ch].max; ch = [1, ch].min; (ch * 255).round }
  end
end

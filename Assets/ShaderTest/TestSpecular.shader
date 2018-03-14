// Upgrade NOTE: replaced '_Object2World' with 'unity_ObjectToWorld'
// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'
// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

// Upgrade NOTE: replaced '_Object2World' with 'unity_ObjectToWorld'

Shader "PBFShaderTest/TestSpecular" {
	Properties{
		_NormalTex("Base (RGB)", 2D) = "white" {}
		_color("Color", Color) = (1, 1, 1, 1)
		_specularColor("Specular Color", Color) = (1, 1, 1, 1) 
		_shininess("Shininess", Int) = 100
		_Cube("Cubemap", Cube) = ""{}
	}
	SubShader{
		GrabPass{ "_GrabTexture" }
		Pass{
			Fog{ Mode off }

			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#pragma target 5.0
			//#pragma enable_d3d11_debug_symbols

			#include "Lighting.cginc"
			#include "UnityCG.cginc"

			float4 _color;
			float4 _specularColor;
			int _shininess;
			samplerCUBE _Cube;
			sampler2D _GrabTexture;

			struct v2f
			{
				float4 pos			: POSITION;
				float2 coord		: TEXCOORD0;
				float3 viewPos		: TEXCOORD1;
				float3 normal		: NORMAL;
				float3 normalDirW 	: TEXCOORD2;
            	float3 viewDirW 	: TEXCOORD3;
			};

			struct fout
			{
				float4 color	: COLOR;
			};


			v2f vert(appdata_base v)
			{
				v2f o;

				o.pos = UnityObjectToClipPos(v.vertex);
				o.coord = v.texcoord.xy;
				//o.modelPos = mul(UNITY_MATRIX_M, v.vertex);
				//o.modelPos = v.vertex;
				o.viewPos = UnityObjectToViewPos(v.vertex);
		
				o.normal = v.normal;
				o.viewDirW = mul(unity_ObjectToWorld, v.vertex).xyz 
               		- _WorldSpaceCameraPos;
            	o.normalDirW = normalize(
               mul(float4(v.normal, 0.0), unity_WorldToObject).xyz);

				return o;
			}

			float4x4 UnityIV()
			{
				//return inverse of unity_matrix_v, which is not provided
				float3x3 subIV = transpose(float3x3(UNITY_MATRIX_V[0].xyz, UNITY_MATRIX_V[1].xyz, UNITY_MATRIX_V[2].xyz));
				//translation
				float3 tIV = (-1 * mul(subIV, float3(UNITY_MATRIX_V[0][3], UNITY_MATRIX_V[1][3], UNITY_MATRIX_V[2][3])));
				return (float4x4(float4(subIV[0], tIV.x), 
								 float4(subIV[1], tIV.y), 
								 float4(subIV[2], tIV.z), 
								 float4(0,0,0,1)));
			}

			fout frag(v2f i)
			{
				fout OUT;

				//float3 eyeSpaceN2 = normalize(mul(UNITY_MATRIX_IT_MV, float4( i.normal, 0.0 ))).xyz;
				float3 eyeSpaceN = normalize( mul( float4( i.normal, 0.0 ), UNITY_MATRIX_IT_MV ).xyz );

				float3 lightDir = normalize(mul(UNITY_MATRIX_V, _WorldSpaceLightPos0)).xyz;
				//float4 lightDir= mul(UNITY_MATRIX_V, _WorldSpaceLightPos0);

				//diffuse
				float diffuseMul = max(0.0, dot(eyeSpaceN, lightDir));
				//float diffuseMul = max(dot(eyeSpaceN, lightDir) * 0.5 + 0.5, 0.0); //wrapped diffuse
				float3 mydiffuse = _color.rgb * diffuseMul * _LightColor0.rgb;

				float3 myambient = _color.rgb * UNITY_LIGHTMODEL_AMBIENT.rgb;  //ShadeSH9(i.normal);

				//row major matrix, inverse of view matrix
				//float3x3 subIV = transpose(float3x3(UNITY_MATRIX_V[0].xyz, UNITY_MATRIX_V[1].xyz, UNITY_MATRIX_V[2].xyz));
				//float3 tV = float3(UNITY_MATRIX_V[0][3], UNITY_MATRIX_V[1][3], UNITY_MATRIX_V[2][3]);
				//float3 tIV = (-1 * mul(subIV, tV));
				//float4x4 IV = float4x4(float4(subIV[0], tIV.x), float4(subIV[1], tIV.y), float4(subIV[2], tIV.z), float4(0,0,0,1));


				float4x4 Id = mul(UnityIV(), UNITY_MATRIX_V);

				//float3x3 id = mul(subVT, subv);

				//specualr 
				float3 viewerVector = normalize(-i.viewPos); //equals to float3(0,0,0) - myPos, zero float3 is camera position

				float3 halfVector = normalize(lightDir + viewerVector);
				float specularMul = pow(max(dot(eyeSpaceN, halfVector), 0.0), _shininess);
				float3 myspecular = _specularColor.rgb * specularMul;

				float3 reflVector = reflect(i.viewDirW, i.normalDirW);
				float3 myreflection = texCUBE(_Cube, reflVector);

				//ShadeSH9 -> unity ambient function, much better than UNITY_LIGHTMODEL_AMBIENT
				//OUT.color = float4(mydiffuse + myambient + myspecular + myreflection, 1.0); //myspecular;

				//float4(abs(UNITY_MATRIX_V[0][3]), abs(UNITY_MATRIX_V[1][3]), abs(UNITY_MATRIX_V[2][3]), 1.0);
				OUT.color = mul(Id, float4(myreflection, 1.0));

				return OUT;
			}
			ENDCG

		}
	}
}

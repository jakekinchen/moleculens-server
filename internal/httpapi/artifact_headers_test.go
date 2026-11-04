package httpapi

import (
	"mime"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"
)

func TestArtifactDownloadsPreserveLegacyHeaders(t *testing.T) {
	for name, wantType := range map[string]string{"Dt.cube": "chemical/x-gaussian-cube", "result.json": "application/json", "structure.sdf": "application/octet-stream"} {
		t.Run(name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), name)
			if err := os.WriteFile(path, []byte("retained scientific bytes"), 0600); err != nil {
				t.Fatal(err)
			}
			response := httptest.NewRecorder()
			serveArtifact(response, httptest.NewRequest("GET", "/download/"+name, nil), path)
			if response.Code != 200 || response.Body.String() != "retained scientific bytes" || response.Header().Get("Content-Type") != wantType {
				t.Fatalf("download contract differs: %d %v %q", response.Code, response.Header(), response.Body.String())
			}
			disposition, params, err := mime.ParseMediaType(response.Header().Get("Content-Disposition"))
			if err != nil || disposition != "attachment" || params["filename"] != name {
				t.Fatalf("download filename/disposition differs: %v %v", response.Header(), err)
			}
		})
	}
}

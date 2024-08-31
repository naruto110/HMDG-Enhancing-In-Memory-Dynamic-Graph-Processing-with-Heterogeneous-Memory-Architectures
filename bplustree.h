#include <time.h>
#include <iostream>
#include <queue>
#include <string>
#include <algorithm>

#include <ext/hash_set>
#include <ext/hash_map>
#define hash_map __gnu_cxx::hash_map
#define hash_set __gnu_cxx::hash_set

#include "avltree.h"

using namespace std;

#define M 1
typedef AVLtree* VAL_TYPE;  

typedef hash_map<int, Tree*> LeafMap;
typedef typename LeafMap::iterator LeafMapIter;


class Inter_Node;

class Node
{
public:
    Node();
    virtual ~Node();
    Node* GetBrother(int& flag);
    Inter_Node* Parent;
    int key[M * 2];//4
    int count;
    int isLeaf;
    void Print();
};
Node::~Node() {}

class Inter_Node : public Node
{
public:
    Inter_Node();
    virtual ~Inter_Node();
    bool Insert(int value, Node* pNode);
    bool Delete(int value);
    int Split(Inter_Node* pNode, int key);
    bool Merge(Inter_Node* pNode);
    bool Slib(Inter_Node* pNode);
    Node* Child[M * 2 + 1];
};
Inter_Node::~Inter_Node()
{
    for (int i = 0; i < M * 2 + 1; i++)
        Child[i] = NULL;
}

class Leaf_Node : public Node
{
public:
    Leaf_Node();
    virtual ~Leaf_Node();
    bool Insert(int value, int neighbor, int tombstone);
    bool Delete(int value);
    int Split(Leaf_Node* pNode);
    bool Merge(Leaf_Node* pNode);
    Leaf_Node* Pre_Node;
    Leaf_Node* Next_Node;

    VAL_TYPE nei_avl[M * 2]; 
    size_t nei_cnt[M * 2] = {0}; 
    // vector<pair<size_t, size_t> > nei_unsorted[M * 2];
};
Leaf_Node::~Leaf_Node() {}

//B+树
class Bplus
{
public:
    Bplus();
    virtual ~Bplus();
    bool Search(int data, string& sPath, VAL_TYPE& nei_avl);
    bool Search_insert(int data, int neighbor, int tombstone);
    bool Insert(int data, int neighbor, int tombstone);  
    bool Delete(int data);
    void Print();
    void treeTraversal();

// protected:
    Leaf_Node* Find(int data);
    bool Add_Node(Inter_Node* pNode, int key, Node* New_Node);
    bool Remove_Node(Inter_Node* pNode, int key);
    Node* Root;
};
Bplus::~Bplus() {}

Node::Node()
{
    isLeaf = true;
    count = 0;
    Parent = NULL;
}

Leaf_Node::Leaf_Node()
{
    isLeaf = true;
    Pre_Node = NULL;
    Next_Node = NULL;

    for (int i=0; i<M * 2; i++)
    {
        nei_avl[i] = NULL;
        nei_cnt[i] = 0;
    }
    // nei_cnt = {0};
}

Inter_Node::Inter_Node()
{
    isLeaf = false;
    for (int i = 0; i < M * 2 + 1; i++)
        Child[i] = NULL;
}

Bplus::Bplus()
{
    Root = NULL;
}

Node* Node::GetBrother(int& flag)
{
    if (NULL == Parent)
        return NULL;

    Node* p = NULL;
    for (int i = 0; i <= Parent->count; i++)
    {
        if (Parent->Child[i] == this)
        {
            if (i == Parent->count) 
            {
                p = Parent->Child[i - 1];
                flag = 1;
            }
            else
            {
                p = Parent->Child[i + 1];
                flag = 2;
            }
        }
    }
    return p;
}

// void Node::Print()
// {
//     for (int i = 0; i < count; i++)
//     {
//         cout << "(" << key[i] << "," << key_words[key[i]] << ")  ";
//         if (i >= count - 1)
//             cout << " | ";
//     }
// }

void Node::Print()
{


}

int Leaf_Node::Split(Leaf_Node* p)
{
    int j = 0;
    for (int i = M; i < M * 2; i++, j++)
    {
        p->key[j] = this->key[i];//this为old node
        p->nei_avl[j] = this->nei_avl[i];
        this->nei_avl[i] = NULL;
        p->nei_cnt[j] = this->nei_cnt[i];
        this->nei_cnt[i] = 0;
    }

    this->count = this->count - j;
    p->count = j;
    return p->key[0];
}

bool Leaf_Node::Delete(int value)
{
    bool found = false;
    int i = 0;
    for (; i < count; i++)
    {
        if (value == key[i])
        {
            found = true;
            break;
        }
    }
    if (false == found)
        return false;
    int j = i;
    for (; j < count - 1; j++)
        key[j] = key[j + 1];
    key[j] = 0;
    count--;
    return true;
}

bool Leaf_Node::Insert(int value, int neighbor, int tombstone)
{
    int i = 0;
    VAL_TYPE tmp_avl = NULL;
    for (; (value > key[i]) && (i < count); i++)
    {
    }
    for (int j = count; j > i; j--)
    {
        key[j] = key[j - 1];
        nei_avl[j] = nei_avl[j - 1];  
        nei_cnt[j] = nei_cnt[j - 1];
    }

    // cout << "test 3: i is " << i << endl;

    key[i] = value;
    tmp_avl = avlInsert(tmp_avl, neighbor, tombstone);

    // cout << "before AVL insert src: " << value << " nei: " << neighbor << " cnt: " << nei_cnt[i] << endl;
    nei_avl[i] = tmp_avl;
    nei_cnt[i] = 1;  

    // cout << "AVL insert src: " << value << " nei: " << neighbor << " cnt: " << nei_cnt[i] << endl;

    // pair<size_t, size_t> tmp_nei(neighbor, tombstone);
    // nei_unsorted[i].push_back(tmp_nei); 

    count++;
    return true;
}

Leaf_Node* Bplus::Find(int data)
{
    int i = 0;
    Node* p = Root;
    Inter_Node* q;  
    while (NULL != p)
    {
        if (p->isLeaf) 
            break;
        for (i = 0; i < p->count; i++) 
        {
            if (data < p->key[i])
                break;
        }
        q = (Inter_Node*)p;
        p = q->Child[i];
    }
    
    return (Leaf_Node*)p;
}

bool Leaf_Node::Merge(Leaf_Node* p)
{
    if (this->count + p->count > M * 2)
        return false;
    for (int i = 0; i < p->count; i++)
        this->Insert(p->key[i], 0, 0);  
    return true;
}

bool Inter_Node::Merge(Inter_Node* p)
{
    key[count] = p->Child[0]->key[0];
    count++;
    Child[count] = p->Child[0];
    for (int i = 0; i < p->count; i++)
    {
        key[count] = p->key[i];
        count++;
        Child[count] = p->Child[i + 1];
    }
    return true;
}

bool Inter_Node::Insert(int value, Node* New)
{
    int i = 0;
    for (; (i < count) && (value > key[i]); i++)
    {
    }
    for (int j = count; j > i; j--)
        key[j] = key[j - 1];
    for (int j = count + 1; j > i + 1; j--)
        Child[j] = Child[j - 1];
    key[i] = value;
    Child[i + 1] = New;
    New->Parent = this;
    count++;
    return true;
}

int Inter_Node::Split(Inter_Node* p, int k)
{
    int i = 0, j = 0;
    if ((k > this->key[M - 1]) && (k < this->key[M]))
    {
        for (i = M; i < M * 2; i++, j++)
            p->key[j] = this->key[i];
        j = 1;
        for (i = M + 1; i <= M * 2; i++, j++)
        {
            this->Child[i]->Parent = p;
            p->Child[j] = this->Child[i];
        }
        this->count = M;
        p->count = M;
        return k;
    }
    int pos = k < this->key[M - 1] ? (M - 1) : M;
    k = this->key[pos];
    j = 0;
    for (i = pos + 1; i < M * 2; i++, j++)
        p->key[j] = this->key[i];
    j = 0;
    for (i = pos + 1; i <= M * 2; i++, j++)
    {
        this->Child[i]->Parent = p;
        p->Child[j] = this->Child[i];
    }
    this->count = pos;
    p->count = M * 2 - pos - 1;
    return k;
}

bool Inter_Node::Delete(int k)
{
    int i = 0;
    for (; (k >= key[i]) && (i < count); i++)
    {
    }
    for (int j = i - 1; j < count - 1; j++)
        key[j] = key[j + 1];
    k = i;
    for (; k < count; k++)
    {
        Child[k] = Child[k + 1];
    }
    Child[k] = NULL;
    count--;
    return true;
}

bool Inter_Node::Slib(Inter_Node* p)
{
    int i, j;
    if (p->key[0] < this->key[0])
    {
        for (i = count; i > 0; i--)
            key[i] = key[i - 1];
        for (j = count + 1; j > 0; j--)
            Child[j] = Child[j - 1];
        key[0] = Child[0]->key[0];
        Child[0] = p->Child[p->count];
    }
    else
    {
        key[count] = p->Child[0]->key[0];
        Child[count + 1] = p->Child[0];
        for (i = 1; i < p->count - 1; i++)
            p->key[i - 1] = p->key[i];
        for (j = 0; j < p->count - 1; j++)
            p->Child[j] = p->Child[j + 1];
    }
    this->count++;
    p->count--;
    return true;
}

bool Bplus::Add_Node(Inter_Node* p, int k, Node* New_Node)  
{
    if (NULL == p || p->isLeaf)
        return false;
    if (p->count < M * 2)
        return p->Insert(k, New_Node);
    Inter_Node* Brother = new Inter_Node;

    int NewKey = p->Split(Brother, k);

   
    if (p->count < Brother->count)
    {
        p->Insert(k, New_Node);
    }
    else if (p->count > Brother->count)
    {
        Brother->Insert(k, New_Node);
    }
    else
    {
        Brother->Child[0] = New_Node;
        New_Node->Parent = Brother;
    }
    Inter_Node* parent = (Inter_Node*)(p->Parent);
    if (NULL == parent)
    {
        parent = new Inter_Node();
        parent->Child[0] = p;
        parent->key[0] = NewKey;
        parent->Child[1] = Brother;
        p->Parent = parent;
        Brother->Parent = parent;
        parent->count = 1;
        Root = parent;
        return true;
    }
    return Add_Node(parent, NewKey, Brother);
}

bool Bplus::Search(int data, string& sPath, VAL_TYPE& tmp_avl)
{
    int i = 0;
    sPath = "The query path: ";
    Node* p = Root;
    if (NULL == p)
        return false;
    Inter_Node* q;
    Leaf_Node* l;
    while (NULL != p)
    {
        if (p->isLeaf)
            break;
        for (i = 0; (i < p->count) && (data >= p->key[i]); i++)
        {
        }
        int k = i > 0 ? i - 1 : i;
        sPath += to_string(p->key[k]);
        sPath += "-->";
        q = (Inter_Node*)p;
        p = q->Child[i];
    }
    if (NULL == p)
        return false;
    sPath += to_string(p->key[0]);
    bool found = false;
    l = (Leaf_Node*) p;  
    for (i = 0; i < l->count; i++)
    {
        if (data == l->key[i])
        {   
            found = true;
            tmp_avl = l->nei_avl[i];
           
            break;
        }
    }
    // cout << "tmp_avl->value: " << tmp_avl->value << endl;
    return found;
}

bool Bplus::Search_insert(int data, int neighbor, int tombstone)
{
    int i = 0;
    Node* p = Root;
    if (NULL == p)
        return false;
    Inter_Node* q;
    Leaf_Node* l;
    while (NULL != p)
    {
        if (p->isLeaf)
            break;
        for (i = 0; (i < p->count) && (data >= p->key[i]); i++)
        {
        }
        int k = i > 0 ? i - 1 : i;
        q = (Inter_Node*)p;
        p = q->Child[i];
    }
    if (NULL == p)
        return false;
    bool found = false;
    l = (Leaf_Node*) p;  
    for (i = 0; i < l->count; i++)
    {
        if (data == l->key[i])
        {   
            found = true;
            
            l->nei_avl[i] = avlInsert(l->nei_avl[i], neighbor, tombstone);
            l->nei_cnt[i] += 1;
            // pair<size_t, size_t> tmp_nei(neighbor, tombstone);
            // l->nei_unsorted[i].push_back(tmp_nei);
            break;
        }
    }
    // cout << "tmp_avl->value: " << tmp_avl->value << endl;
    return found;
}

bool Bplus::Insert(int data, int neighbor, int tombstone) 
{
    
    string a;
    VAL_TYPE tmp_avl = NULL;
    if (true == Search_insert(data, neighbor, tombstone))
    {
        return true;      
    }
    
    Leaf_Node* Old_Node = Find(data);
   
    if (NULL == Old_Node) 
    {
        Old_Node = new Leaf_Node;
        Root = Old_Node;
    }
    

    if (Old_Node->count < M * 2) 
    {
        // cout << "test 2: insert the first key: " << data << endl;
        return Old_Node->Insert(data, neighbor, tombstone);  
    }

    Leaf_Node* New_Node = new Leaf_Node;
    
    int k = Old_Node->Split(New_Node);

    Leaf_Node* OldNext = Old_Node->Next_Node;  
    Old_Node->Next_Node = New_Node;
    New_Node->Next_Node = OldNext;
    New_Node->Pre_Node = Old_Node;

    if (NULL != OldNext)
        OldNext->Pre_Node = New_Node;

    if (data < k)
    {
        Old_Node->Insert(data, neighbor, tombstone);
    }
    else
    {
        New_Node->Insert(data, neighbor, tombstone);
    }
    Inter_Node* parent = (Inter_Node*)(Old_Node->Parent);

    if (NULL == parent)
    {
        Inter_Node* New_Root = new Inter_Node;
        New_Root->Child[0] = Old_Node;
        New_Root->key[0] = k;
        New_Root->Child[1] = New_Node;
        Old_Node->Parent = New_Root;
        New_Node->Parent = New_Root;
        New_Root->count = 1;
        Root = New_Root;
        return true;
    }

    return Add_Node(parent, k, New_Node);
}

bool Bplus::Delete(int data)
{
    Leaf_Node* Old_Node = Find(data); 
    if (NULL == Old_Node)
        return false;
    if (false == Old_Node->Delete(data)) 
        return false;
    Inter_Node* parent = (Inter_Node*)(Old_Node->Parent);
    if (NULL == parent)
    {
        if (0 == Old_Node->count)
        {
            delete Old_Node;
            Root = NULL;
        }
        return true;
    }
    if (Old_Node->count >= M)
    {
        for (int i = 0; (i < parent->count) && (data >= parent->key[i]); i++)
        {
            if (parent->key[i] == data)
                parent->key[i] = Old_Node->key[0];
        }
        return true;
    }
    
    int flag = 1;
    Leaf_Node* Brother = (Leaf_Node*)(Old_Node->GetBrother(flag));
    int NewData = 0;
    if (Brother->count > M)
    {
        if (1 == flag)
        {
            NewData = Brother->key[Brother->count - 1];
        }
        else
        {
            NewData = Brother->key[0];
        }
        Old_Node->Insert(NewData, 0, 0);
        Brother->Delete(NewData);
        
        if (1 == flag)
        {
            for (int i = 0; i <= parent->count; i++)
            {
                if (parent->Child[i] == Old_Node && i > 0)
                    parent->key[i - 1] = Old_Node->key[0];
            }
        }
        else
        {
            for (int i = 0; i <= parent->count; i++)
            {
                if (parent->Child[i] == Old_Node && i > 0)
                    parent->key[i - 1] = Old_Node->key[0];
                if (parent->Child[i] == Brother && i > 1)
                    parent->key[i - 1] = Brother->key[0];
            }
        }
        return true;
    }
    int NewKey = 0;
    if (1 == flag)
    {
        Brother->Merge(Old_Node);
        NewKey = Old_Node->key[0];
        Leaf_Node* OldNext = Old_Node->Next_Node;
        Brother->Next_Node = OldNext;
        if (NULL != OldNext)
            OldNext->Pre_Node = Brother;
        delete Old_Node;
    }
    else
    {
        Old_Node->Merge(Brother);
        NewKey = Brother->key[0];
        Leaf_Node* OldNext = Brother->Next_Node;
        Old_Node->Next_Node = OldNext;
        if (NULL != OldNext)
            OldNext->Pre_Node = Old_Node;
        delete Brother;
    }
    return Remove_Node(parent, NewKey);
}

bool Bplus::Remove_Node(Inter_Node* p, int k)
{
    if (false == p->Delete(k))
    {
        return false;
    }
    Inter_Node* parent = (Inter_Node*)(p->Parent);
    if (NULL == parent)
    {
        if (0 == p->count)
        {
            Root = p->Child[0];
            delete p;
        }
        return true;
    }
    if (p->count >= M)
    {
        
        for (int i = 0; (i < parent->count) && (k >= parent->key[i]); i++)
        {
            if (parent->key[i] == k)
                parent->key[i] = p->key[0];
        }
        return true;
    }
    
    int flag = 1;
    Inter_Node* Brother = (Inter_Node*)(p->GetBrother(flag));
    if (Brother->count > M)
    {
        p->Slib(Brother);
        if (1 == flag)
        {
            for (int i = 0; i <= parent->count; i++)
            {
                if (parent->Child[i] == p && i > 0)
                    parent->key[i - 1] = p->key[0];
            }
        }
        else
        {
            for (int i = 0; i <= parent->count; i++)
            {
                if (parent->Child[i] == p && i > 0)
                    parent->key[i - 1] = p->key[0];
                if (parent->Child[i] == Brother && i > 0)
                    parent->key[i - 1] = Brother->key[0];
            }
        }
        return true;
    }
    
    int NewKey = 0;
    if (1 == flag)
    {
        Brother->Merge(p);
        NewKey = p->key[0];
        delete p;
    }
    else
    {
        p->Merge(Brother);
        NewKey = Brother->key[0];
        delete Brother;
    }
    return Remove_Node(parent, NewKey);
}

void Bplus::Print()
{
    Node* p = Root;
    if (NULL == p)
        return;
    Inter_Node* a;
    int H = 0;
    queue<Node*> q;
    queue<int> h;
    q.push(p);
    h.push(1);
    while (!q.empty())
    {
        p = q.front();
        if (H != h.front())
        {
            cout << endl;
            cout << H << endl;
            H = h.front();
        }
        q.pop();
        h.pop();
        p->Print();
        if (NULL != p && !p->isLeaf)
        {
            a = (Inter_Node*)p;
            for (int i = 0; i <= p->count; i++)
            {
                q.push(a->Child[i]);
                h.push(H + 1);
            }
        }
    }
}

void Bplus::treeTraversal()
{
    Node* p = Root;
    if (NULL == p)
        return;
    
    Inter_Node* q;
    Leaf_Node* l;
    while (NULL != p) 
    {
        if (p->isLeaf)
            break;
        q = (Inter_Node*)p;
        p = q->Child[0];  
    }
    l = (Leaf_Node*) p; 
    while (l != NULL)
    {
        for (size_t keyIdx=0; keyIdx<l->count; keyIdx++)
        {
            size_t src = l->key[keyIdx];
            VAL_TYPE tmp_avl;
            tmp_avl = l->nei_avl[keyIdx];
            cout << "src: " << src << " neighbor list in DRAM: ";
            displayTree_Mid(tmp_avl);
            cout << endl;
        }
        l = l->Next_Node;
    }
    
}

void test_generate(Bplus* bplus)
{
    int a[] = {1, 1, 1, 1, 1, 1, 2, 3};
    int b[] = {2, 3, 6, 7, 8, 9, 4, 5};
    int c[] = {1, 1, 1, 1, 1, 1, 1, 0};

    for (int i = 0; i < 8; i++)
    {
        // cout << "insert " << i << " begin." << endl;
        bplus->Insert(a[i], b[i], c[i]);
        // cout << "insert " << i << " over." << endl;
    }
    
    cout << "Create B+ tree suceess!" << endl;
    cout << endl;
}

void query_neighbor(Bplus* bplus, int src)
{
    string sPath;
    VAL_TYPE tmp_avl;
    bplus->Search(src, sPath, tmp_avl);
    cout << "src: " << src << " neighbor list in DRAM: ";
    displayTree_Mid(tmp_avl);
    cout << endl;
}

void insert_edge(Bplus* bplus, int src, int des, int tombstone)
{
    bool success = bplus->Insert(src, des, tombstone);
    // if (true == success)
    // {
    //     cout << "OK" << endl;
    // }
    // else
    // {
    //     cout << "FAIL" << endl;
    // }
    // cout << endl;
}

void test4_genavl(void)
{
    VAL_TYPE tmp_avl = NULL;
    int a[] = {1, 1, 1, 1, 1, 1, 2, 3};
    int b[] = {2, 3, 6, 7, 8, 9, 4, 5};
    int c[] = {1, 1, 1, 1, 1, 1, 1, 0};

    for (int i = 0; i < 8; i++)
    {
        tmp_avl = avlInsert(tmp_avl, b[i], c[1]);
        // cout << "tmp_val value: " << tmp_avl->value << " b[i]: " << b[i] << endl;
    }
    cout << "tmp_val value: " << tmp_avl->value << " depth: " << tmp_avl->depth << endl;
    displayTree_Mid(tmp_avl);
    cout << endl;
}

void test5(Bplus* bplus)
{
    bplus->Print();
    cout << endl;
}

